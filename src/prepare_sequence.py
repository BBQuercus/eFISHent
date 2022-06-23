"""
Download Entrez Gene Probes from NCBI.
Select the right strand and select intronic/exonic regions.
"""

import logging
import os
import subprocess

import Bio.SeqIO
import luigi

from config import GeneralConfig
from config import SequenceConfig
import util


class DownloadEntrezGeneSequence(luigi.Task):
    """Download genomic sequence as fasta from NCBI."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        fname = f"{util.get_gene_name()}_entrez.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def run(self):
        has_ensembl = SequenceConfig().ensemble_id
        has_gene_and_organism = (
            SequenceConfig().gene_name and SequenceConfig().organism_name
        )
        if not has_ensembl and not has_gene_and_organism:
            raise ValueError(
                "For downloading Entrez Gene Probes, "
                " you need to specify the gene name and organism name or provide an Emsembl ID."
            )

        if has_ensembl:
            query = SequenceConfig().ensemble_id
        else:
            query = f"{SequenceConfig().gene_name} [GENE] {SequenceConfig().organism_name} [ORGN]"
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
        args_entrez = [*args_search, "|", *args_link, "|", *args_fetch]
        self.logger.debug(f"Fetching from Entrez using query '{query}'.")
        fasta = subprocess.check_output(args_entrez, stderr=subprocess.STDOUT).decode()

        if not fasta:
            raise LookupError(
                "No Entrez gene found. "
                "Please check your gene name and organism and try again. "
                "Or provide a custom fasta file as `sequence-file`."
            )

        with self.output().open("w") as outfile:
            outfile.write(fasta)


class BuildBlastDatabase(luigi.Task):
    def output(self):
        return [
            luigi.LocalTarget(f"{util.get_genome_name()}.{extension}")
            for extension in ["ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"]
        ]

    def run(self):
        args_blast = [
            "makeblastdb",
            "-dbtype",
            "nucl",
            "-in",
            GeneralConfig().reference_genome,
            "-out",
            util.get_genome_name(),
        ]
        subprocess.check_call(args_blast)


class PrepareSequence(luigi.Task):
    """Prepare gene sequence.

    Select right strand and exon/intron/both."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        tasks = {}
        if SequenceConfig().sequence_file == "None":
            tasks["entrez"] = DownloadEntrezGeneSequence()
        if SequenceConfig().is_exonic or SequenceConfig().is_intronic:
            tasks["blast"] = BuildBlastDatabase()
        return tasks

    def output(self):
        fname = f"{util.get_gene_name()}_sequence.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def select_gene_region(
        self, sequence: Bio.SeqRecord.SeqRecord
    ) -> Bio.SeqRecord.SeqRecord:
        """Select exon/intronic regions."""
        # TODO select exon/intron/both
        # Create blast database
        # Blast sequence against database
        # Parse blast output
        # Select regions based on genome annotation

    def run(self):
        input_file = (
            self.input()["entrez"].path
            if "entrez" in self.input()
            else SequenceConfig().sequence_file
        )
        sequence = list(Bio.SeqIO.parse(input_file, format="fasta"))

        if not sequence:
            raise ValueError(
                "No records found in fasta file. "
                "Please ensure at least one sequence is present."
            )
        if len(sequence) > 1:
            self.logger.warning(
                f"{util.UniCode.warn} More than one record (potential isoforms) found in fasta file. "
                "Will default to take the first sequence."
            )
        sequence = sequence[0]

        if SequenceConfig().is_plus_strand:
            sequence = sequence.reverse_complement(id=True, description=True)
            self.logger.debug("Converted sequence to reverse complement.")

        # sequence = self.select_gene_region(sequence)

        Bio.SeqIO.write(sequence, self.output().path, format="fasta")
