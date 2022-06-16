"""
Create bowtie index for a given genome.
Align probes to the reference genome.
Filter probes based on alignment score and uniqueness.
"""

import logging
import os

from tqdm import tqdm
import Bio.SeqIO
import gtfparse
import luigi
import numpy as np
import pandas as pd
import pysam

from config import GeneralConfig
from config import ProbeConfig
from config import SequenceConfig
from basic_filtering import BasicFiltering
from prepare_sequence import BuildBlastDatabase
from prepare_sequence import PrepareSequence
import constants
import util


class BuildBowtieIndex(luigi.Task):
    """Create bowtie index for a given reference genome."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return [
            luigi.LocalTarget(util.get_genome_name() + ext)
            for ext in [
                ".1.ebwt",
                ".2.ebwt",
                ".3.ebwt",
                ".4.ebwt",
                ".rev.1.ebwt",
                ".rev.2.ebwt",
            ]
        ]

    def run(self):
        os.system(
            f"bowtie-build {os.path.abspath(GeneralConfig().reference_genome)} {util.get_genome_name()} "
            f"--threads {GeneralConfig().threads}"
        )


class AlignProbeCandidates(luigi.Task):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")
    is_blast_required = (
        ProbeConfig().encode_count_table is not None and SequenceConfig().is_endogenous
    )

    def _requires(self):
        tasks = [BuildBowtieIndex()]
        if self.is_blast_required:
            tasks.append(BuildBlastDatabase())
        return luigi.task.flatten(
            [
                *tasks,
                super(AlignProbeCandidates, self)._requires(),
            ]
        )

    def requires(self):
        tasks = {"probes": BasicFiltering()}
        # if self.is_blast_required:
        #     tasks.append(PrepareSequence())
        return tasks

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_aligned.fasta")
        )

    @staticmethod
    def align_probes(fname_fasta: str, fname_sam: str) -> None:
        """Align probes to the reference genome."""
        # Convert fasta to fastq - bowtie doesn't return read names if not fastq...
        fname_fastq = fname_fasta.rstrip("a") + "q"
        os.system(f"seqtk seq -F 'I' {fname_fasta} > {fname_fastq}")

        # TODO change params to match endo/non-endo (currently only endo!)
        # TODO --quiet output to logger
        endo_exo_param = (
            ""
            if SequenceConfig().is_endogenous
            else f"-m {ProbeConfig().max_off_targets}"
        )
        os.system(
            f"bowtie\
                -x {util.get_genome_name()} {fname_fastq}\
                --threads {GeneralConfig().threads}\
                {endo_exo_param}\
                --sam -S {fname_sam}"
        )

    @staticmethod
    def filter_unique_probes(fname_sam: str) -> pd.DataFrame:
        """Filter probes based on alignment score and uniqueness."""
        # 60 - uniquely mapped read, regardless of number of mismatches / indels
        # 4 – flag for unmapped read
        pysam_parameters = (
            ["--min-MQ", "60"]
            if SequenceConfig().is_endogenous
            else ["--require_flags", "4"]
        )
        filtered_sam = pysam.view(*pysam_parameters, fname_sam)

        # Parse tab and newline delimited pysam output
        df = pd.DataFrame(
            [row.split("\t") for row in filtered_sam.split("\n")],
            columns=constants.SAMFILE_COLUMNS,
        )
        # order = list(df[df["qname"] != ""].groupby("qname").size().sort_values().keys())
        # sequences_ordered = []
        # for name in order:
        #     for seq in sequences:
        #         if seq.id == name:
        #             sequences_ordered.append(seq)
        return df

    def filter_gene_of_interest(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter probes based on gene of interest."""
        fname = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_blast.fasta"
        )
        os.system(
            f"blastn\
                -db {util.get_genome_name()}\
                -query {self.input()[1].path}\
                -task megablast\
                -outfmt 6\
                -num_threads {GeneralConfig().threads} > {fname}"
        )
        df_blast = pd.read_csv("./test.txt", sep="\t", names=constants.BLAST_COLUMNS)
        df_blast = df_blast[df_blast["evalue"] == 0]
        # TODO add filter for gene of interest
        return df

    def filter_fpkm(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter candidates based on FPKM."""
        # Read gtf file
        df_gtf = gtfparse.read_gtf(GeneralConfig().reference_annotation)

        # Read count table
        df_counts = pd.read_csv(ProbeConfig().encode_count_table, sep="\t")
        counts = df_counts[df_counts["gene_id"].str.contains("ENSG")].copy()
        counts["clean_id"] = counts["gene_id"].apply(lambda x: x.split(".")[0])

        # Merge gtf and count table to get FPKM at each genomic location
        df_fpkm = pd.merge(
            counts[constants.COUNTS_COLUMNS],
            df_gtf[constants.GTF_COLUMNS],
            how="left",
            left_on="clean_id",
            right_on="gene_id",
        ).dropna()

        if self.is_blast_required:
            df_fpkm = self.filter_gene_of_interest(df_fpkm)

        # Match FPKM to the probe candidates' genomic locations
        dfs = []
        for sequence, df_sequence in tqdm(df.groupby("rname")):
            low = df_fpkm[df_fpkm["seqname"] == sequence]["start"].astype(int).values
            high = df_fpkm[df_fpkm["seqname"] == sequence]["end"].astype(int).values
            vals = df_sequence["pos"].astype(int).values
            i, j = np.where((vals[:, None] >= low) & (vals[:, None] <= high))
            dfs.append(
                pd.concat(
                    [
                        df_sequence.reset_index(drop=True)
                        .iloc[i]
                        .reset_index(drop=True),
                        df_fpkm[df_fpkm["seqname"] == sequence]
                        .reset_index(drop=True)
                        .iloc[j]
                        .reset_index(drop=True),
                    ],
                    axis=1,
                )
            )
        dfs = pd.concat(dfs, ignore_index=True)

        # Filter out probes with FPKM > max_fpkm
        df_max_fpkm = dfs.groupby("qname", as_index=False)["FPKM"].max()
        candidates = df_max_fpkm[
            df_max_fpkm["FPKM"] > ProbeConfig().max_off_target_fpkm
        ]
        return candidates

    def run(self):
        fname_fasta = self.input()["probes"].path
        fname_sam = os.path.splitext(fname_fasta)[0] + ".sam"
        self.align_probes(fname_fasta, fname_sam)
        df = self.filter_unique_probes(fname_sam)

        if ProbeConfig().encode_count_table != "None":
            df = self.filter_fpkm(df)

        # Filter candidates found in samfile
        sequences = list(Bio.SeqIO.parse(fname_fasta, "fasta"))
        candidates = [seq for seq in tqdm(sequences) if seq.id in df["qname"].values]
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, "fasta")
