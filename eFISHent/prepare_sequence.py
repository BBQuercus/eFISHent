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

    def get_entrez_query(
        self, ensembl_id: str, gene_name: str, organism_name: str
    ) -> str:
        """Retrieve the query depending on sequence configuration arguments passed."""
        if ensembl_id:
            return f"({ensembl_id})" + (
                f" AND {organism_name}[Organism]" if organism_name else ""
            )

        if gene_name and organism_name:
            return f"({gene_name}[Gene Name]) AND {organism_name}[Organism]"

        raise ValueError(
            "For downloading Entrez Gene Probes, "
            " you need to specify the gene name and organism name or provide an Emsembl ID."
        )

    def fetch_entrez(self, query: str) -> str:
        """Retrieve the fasta sequence data from entrez."""
        args_search = ["esearch", "-db", "gene", "-query", query]
        args_link = [
            "elink",
            "-db",
            "gene",
            "-target",
            "nuccore",
            "-name",
            "gene_nuccore_refseqrna",
        ]
        args_fetch = ["efetch", "-format", "fasta"]
        self.logger.debug(f"Fetching from Entrez using query '{query}'.")

        search = subprocess.run(args_search, check=True, capture_output=True)
        link = subprocess.run(args_link, input=search.stdout, capture_output=True)
        fetch = subprocess.run(args_fetch, input=link.stdout, capture_output=True)
        fasta = fetch.stdout.decode()

        if not fasta:
            raise LookupError(
                "No Entrez gene found. "
                "Please check your gene name and organism and try again. "
                "Or provide a custom fasta file as `sequence-file`."
            )
        return fasta

    def run(self):
        if not SequenceConfig().is_plus_strand:
            self.logger.warning(
                (
                    "Downloading from entrez but using negative strand as target! "
                    "Entrez usually provides sequences 5'-3'. "
                    "Proceed with caution."
                )
            )
        query = self.get_entrez_query(
            SequenceConfig().ensembl_id,
            SequenceConfig().gene_name,
            SequenceConfig().organism_name,
        )
        fasta = self.fetch_entrez(query)

        with self.output().open("w") as outfile:
            outfile.write(fasta)
        self.logger.info(f"Downloaded sequence from entrez using query '{query}'.")


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
        """Select a single sequence as template."""
        if not sequences:
            raise ValueError(
                "No records found in fasta file. "
                "Please ensure at least one sequence is present."
            )
        if len(sequences) > 1:
            self.logger.warning(
                f"{util.UniCode.warn} More than one record (potential isoforms) found in fasta file. "
                "Will default to take the first sequence."
            )
        return sequences[0]

    def select_strand(
        self, sequence: Bio.SeqRecord.SeqRecord, is_plus_strand: bool
    ) -> Bio.SeqRecord.SeqRecord:
        """Select right strand for probes to target."""
        if is_plus_strand:
            sequence = sequence.reverse_complement(id=True, description=True)
            self.logger.debug("Converted sequence to reverse complement.")
        return sequence

    def run(self):
        input_file = (
            self.input()["entrez"].path
            if "entrez" in self.input()
            else SequenceConfig().sequence_file
        )
        sequences = list(Bio.SeqIO.parse(input_file, format="fasta"))
        sequence = self.select_sequence(sequences)
        sequence = self.select_strand(sequence, SequenceConfig().is_plus_strand)
        self.logger.info(f'Selected sequence and strand called "{sequence.id}".')
        Bio.SeqIO.write(sequence, self.output().path, format="fasta")
