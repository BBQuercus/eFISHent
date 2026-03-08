"""Download Entrez Gene Probes from NCBI.

Select the right strand and select intronic/exonic regions.
"""

from typing import List
import logging
import os
import subprocess

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .config import SequenceConfig


class DownloadEntrezGeneSequence(luigi.Task):
    """Download genomic sequence as fasta from NCBI."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        fname = f"{util.get_gene_name()}_entrez.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def get_nuccore_query(
        self, ensembl_id: str, gene_name: str, organism_name: str
    ) -> str:
        """Build an NCBI nuccore query that returns only RefSeq transcripts.

        Searches the nuccore database directly with RefSeq and biomol filters
        to avoid returning genomic sequences, ESTs, or other non-transcript records.
        For Ensembl IDs, first resolves to an NCBI Gene ID via the gene database.
        """
        refseq_filter = "refseq[filter] AND (biomol_mrna[prop] OR biomol_rna[prop])"

        if ensembl_id:
            # Resolve Ensembl ID to NCBI Gene ID via the gene database
            ensembl_query = f"({ensembl_id})"
            if organism_name:
                ensembl_query += f" AND {organism_name}[Organism]"

            search = subprocess.run(
                ["esearch", "-db", "gene", "-query", ensembl_query],
                check=True,
                capture_output=True,
            )
            uid_result = subprocess.run(
                ["efetch", "-format", "uid"],
                input=search.stdout,
                check=True,
                capture_output=True,
            )
            gene_id = uid_result.stdout.decode().strip()
            if not gene_id:
                raise LookupError(
                    f"No NCBI Gene found for Ensembl ID '{ensembl_id}'. "
                    "Please check the ID and try again."
                )
            # Take first Gene ID if multiple returned
            gene_id = gene_id.split("\n")[0].strip()
            self.logger.debug(f"Resolved Ensembl ID '{ensembl_id}' to Gene ID {gene_id}.")
            return f"{gene_id}[Gene ID] AND {refseq_filter}"

        if gene_name and organism_name:
            return (
                f"({gene_name}[Gene Name]) AND {organism_name}[Organism]"
                f" AND {refseq_filter}"
            )

        raise ValueError(
            "For downloading Entrez Gene Probes, "
            "you need to specify the gene name and organism name or provide an Ensembl ID."
        )

    def fetch_entrez(
        self, ensembl_id: str, gene_name: str, organism_name: str
    ) -> str:
        """Retrieve the fasta sequence data from NCBI.

        Queries the nuccore database directly with RefSeq transcript filters
        to ensure only curated mRNA/ncRNA sequences are returned.
        """
        query = self.get_nuccore_query(ensembl_id, gene_name, organism_name)
        self.logger.debug(f"Fetching from Entrez nuccore using query '{query}'.")

        search = subprocess.run(
            ["esearch", "-db", "nuccore", "-query", query],
            check=True,
            capture_output=True,
        )
        fetch = subprocess.run(
            ["efetch", "-format", "fasta"],
            input=search.stdout,
            check=True,
            capture_output=True,
        )
        fasta = fetch.stdout.decode()

        # Check for empty results or missing FASTA content
        if not fasta or ">" not in fasta:
            raise LookupError(
                "No RefSeq transcript found for your query. "
                "Please check your gene name and organism and try again. "
                "Or provide a custom fasta file with --sequence-file."
            )
        return fasta

    def run(self):
        util.log_stage_start(self.logger, "DownloadEntrezGeneSequence")
        if not SequenceConfig().is_plus_strand:
            self.logger.warning(
                (
                    "Downloading from entrez but using negative strand as target! "
                    "Entrez usually provides sequences 5'-3'. "
                    "Proceed with caution."
                )
            )
        fasta = self.fetch_entrez(
            SequenceConfig().ensembl_id,
            SequenceConfig().gene_name,
            SequenceConfig().organism_name,
        )

        # Log which sequences were retrieved
        headers = [line for line in fasta.split("\n") if line.startswith(">")]
        n_records = len(headers)
        if n_records == 1:
            self.logger.info(f"Downloaded 1 transcript: {headers[0][1:]}")
        else:
            self.logger.info(f"Downloaded {n_records} transcript isoforms from NCBI:")
            for header in headers:
                self.logger.info(f"  {header[1:]}")

        with self.output().open("w") as outfile:
            outfile.write(fasta)


class PrepareSequence(luigi.Task):
    """Prepare gene sequence.

    Select right strand and exon/intron/both.
    """

    logger = logging.getLogger("custom-logger")

    def requires(self):
        tasks = {}
        if not SequenceConfig().sequence_file:
            tasks["entrez"] = DownloadEntrezGeneSequence()
        return tasks

    def output(self):
        fname = f"{util.get_gene_name()}_sequence.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def select_sequence(
        self, sequences: List[Bio.SeqRecord.SeqRecord]
    ) -> Bio.SeqRecord.SeqRecord:
        """Select a single sequence as template.

        When multiple records are present (e.g. transcript isoforms),
        selects the longest sequence to maximize probe coverage.
        """
        if not sequences:
            raise ValueError(
                "No records found in fasta file. "
                "Please ensure at least one sequence is present."
            )
        if len(sequences) > 1:
            longest = max(sequences, key=lambda s: len(s.seq))
            self.logger.info(
                f"Found {len(sequences)} transcript isoforms. "
                f"Selected longest: {longest.id} ({len(longest.seq):,} nt)."
            )
            return longest
        return sequences[0]

    def select_strand(
        self, sequence: Bio.SeqRecord.SeqRecord, is_plus_strand: bool
    ) -> Bio.SeqRecord.SeqRecord:
        """Select right strand for probes to target."""
        if is_plus_strand:
            sequence = sequence.reverse_complement(id=True, description=True)
            self.logger.debug("Converted sequence to reverse complement.")
        return sequence

    @staticmethod
    def is_fasta_formatted(fname):
        if not os.path.exists(fname):
            return False

        with open(fname, "r") as fasta_file:
            in_sequence = False
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    in_sequence = True
                elif not in_sequence or not all(
                    base in "ATCGNatcgn" for base in line
                ):
                    return False
            return True

    def run(self):
        util.log_stage_start(self.logger, "PrepareSequence")
        input_file = (
            self.input()["entrez"].path
            if "entrez" in self.input()
            else SequenceConfig().sequence_file
        )
        if not self.is_fasta_formatted(input_file):
            raise ValueError(
                "Fasta file is incorrectly formatted. "
                "Check wikipedia: https://en.wikipedia.org/wiki/FASTA_format"
            )
        sequences = list(Bio.SeqIO.parse(input_file, format="fasta"))
        sequence = self.select_sequence(sequences)
        sequence = self.select_strand(sequence, SequenceConfig().is_plus_strand)
        self.logger.info(f'Selected sequence and strand called "{sequence.id}".')
        Bio.SeqIO.write(sequence, self.output().path, format="fasta")
