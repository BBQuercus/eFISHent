"""
Download Entrez Gene Probes from NCBI.
Select the right strand and select intronic/exonic regions.
"""

import os
import logging

import luigi
import Bio.SeqIO

from config import GeneralConfig
from config import SequenceConfig
import util


class DownloadEntrezGeneSequence(luigi.Task):
    """Download genomic sequence as fasta from NCBI."""

    logger = logging.getLogger("luigi-interface")

    def output(self):
        fname = f"{util.get_gene_name()}_entrez.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def run(self):
        if not SequenceConfig().gene_name or not SequenceConfig().organism_name:
            raise ValueError(
                "For downloading Entrez Gene Probes, "
                " you need to specify the gene name and organism name."
            )

        fasta = os.popen(
            f"esearch\
                -db gene\
                -query '{SequenceConfig().gene_name} [GENE]\
                    {SequenceConfig().organism_name} [ORGN]' |"
            " elink\
                -db gene\
                -target nuccore\
                -name gene_nuccore_refseqrna |"
            " efetch -format fasta"
        ).read()
        self.logger.debug(
            f"Fetched Entrez Gene Probes for {SequenceConfig().gene_name} in {SequenceConfig().organism_name}."
        )

        if not fasta:
            raise LookupError(
                "No Entrez Gene Probes found. "
                "Please check your gene name and organism and try again."
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
        os.system(
            f"makeblastdb\
                -dbtype nucl\
                -in {GeneralConfig().reference_genome}\
                -out {util.get_genome_name()}"
        )


class PrepareSequence(luigi.Task):
    """Prepare gene sequence.

    Select right strand and exon/intron/both."""

    logger = logging.getLogger("luigi-interface")

    def requires(self):
        tasks = {}
        if not SequenceConfig().sequence_file:
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
            if self.input()["entrez"]
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
                "More than one record (potential isoforms) found in fasta file. "
                "Will default to take the first sequence."
            )
            sequence = sequence[0]

        if SequenceConfig().is_plus_strand:
            sequence = sequence.reverse_complement(id=True, description=True)
            self.logger.debug("Converted sequence to reverse complement.")

        # sequence = self.select_gene_region(sequence)

        Bio.SeqIO.write(sequence, self.output().path, format="fasta")
