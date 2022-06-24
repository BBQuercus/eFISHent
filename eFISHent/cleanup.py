"""
Clean up the output directory and remove the files that are not needed.
Prettify the output files.

Output files:
    - .fasta with the probe sequences
    - .csv with probe information (Tm, GC, deltaG, etc.)
"""

import glob
import logging
import os

import luigi
import pandas as pd
import Bio.SeqIO

from .alignment import AlignProbeCandidates
from .basic_filtering import get_gc_content
from .basic_filtering import get_melting_temp
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .kmers import BuildJellyfishIndex
from .kmers import get_max_kmer_count
from .optimization import OptimizeProbeCoverage
from .secondary_structure import get_free_energy
from . import util


class CleanUpOutput(luigi.Task):
    """Clean up the output files and remove the intermediaries that are not needed."""

    logger = logging.getLogger("custom-logger")
    is_blast_required = (
        ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
    )

    def requires(self):
        tasks = {
            "optimize": OptimizeProbeCoverage(),
            "jellyfish": BuildJellyfishIndex(),
        }
        if self.is_blast_required:
            tasks["alignment"] = AlignProbeCandidates()
        return tasks

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("fasta", f"{util.get_gene_name()}.fasta"),
                ("table", f"{util.get_gene_name()}.csv"),
            ]
        }

    def prettify_table(self, sequences: list):
        """Create table with probe information."""
        # Create basic table with probe start/end positions and name
        df = util.create_data_table(sequences)

        # Add data columns
        df["GC"] = [round(get_gc_content(seq.seq), 2) for seq in sequences]
        df["TM"] = [
            round(
                get_melting_temp(
                    seq.seq,
                    ProbeConfig().na_concentration,
                    ProbeConfig().formamide_concentration,
                ),
                2,
            )
            for seq in sequences
        ]
        df["deltaG"] = [get_free_energy(seq) for seq in sequences]

        # Add jellyfish count
        df["kmers"] = [
            get_max_kmer_count(seq, self.input()["jellyfish"].path) for seq in sequences
        ]

        # Add max FPKM
        if self.is_blast_required:
            df_fpkm = pd.read_csv(self.input()["alignment"]["table"].path)
            df["FPKM"] = pd.merge(
                df, df_fpkm, how="left", left_on="name", right_on="qname"
            )["FPKM"].fillna(0)

        df["name"] = [
            f"{util.get_gene_name(hashed=False)}-{idx + 1}" for idx in df.index
        ]
        return df

    def prettify_sequences(self, df: pd.DataFrame) -> list:
        """Clean up sequence id's and descriptions."""
        return [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(seq), id=name, name=name, description=""
            )
            for seq, name in zip(df["sequence"], df["name"])
        ]

    def run(self):
        intermediary_files = glob.glob(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}*")
        )

        sequences = list(
            Bio.SeqIO.parse(self.input()["optimize"]["probes"].path, "fasta")
        )
        df = self.prettify_table(sequences)
        sequences = self.prettify_sequences(df)

        df.to_csv(self.output()["table"].path, index=False)
        util.log_and_check_candidates(self.logger, "CleanUpOutput", len(sequences))
        Bio.SeqIO.write(sequences, self.output()["fasta"].path, "fasta")

        if not RunConfig().save_intermediates:
            for filename in intermediary_files:
                if filename.endswith(
                    (".fasta", ".sam", ".fastq", ".png", ".txt", ".csv")
                ):
                    self.logger.debug(f"Removing {filename}")
                    os.remove(filename)
