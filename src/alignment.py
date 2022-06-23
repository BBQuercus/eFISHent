"""
Create bowtie index for a given genome.
Align probes to the reference genome.
Filter probes based on alignment score and uniqueness.
"""

import logging
import os
import subprocess
import warnings

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
        args_bowtie = [
            "bowtie-build",
            os.path.abspath(GeneralConfig().reference_genome),
            util.get_genome_name(),
            "--threads",
            str(GeneralConfig().threads),
        ]
        subprocess.check_call(args_bowtie)


class PrepareAnnotationFile(luigi.Task):
    """Convert gtf annotation file to parquet.

    Doesn't change any data but greatly speeds up subsequent I/O (~10x increase).
    """

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(GeneralConfig().reference_annotation + ".parq")

    def run(self):
        warnings.filterwarnings("ignore")
        df = gtfparse.read_gtf(GeneralConfig().reference_annotation)
        warnings.filterwarnings("default")
        self.logger.debug("Read gtf file with gtfparse")
        df.to_parquet(self.output().path)


class AlignProbeCandidates(luigi.Task):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")
    is_blast_required = (
        ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
    )

    def _requires(self):
        """Hidden / not directly accessed requirements."""
        tasks = [BuildBowtieIndex()]
        if self.is_blast_required:
            tasks.append(BuildBlastDatabase())
        return luigi.task.flatten(
            [*tasks, super(AlignProbeCandidates, self)._requires()]
        )

    def requires(self):
        tasks = {"probes": BasicFiltering()}
        if self.is_blast_required:
            tasks["blastseq"] = PrepareSequence()
            tasks["gtf"] = PrepareAnnotationFile()
        return tasks

    def output(self):
        tasks = {
            "fasta": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_aligned.fasta"
                )
            )
        }
        if self.is_blast_required:
            tasks["table"] = luigi.LocalTarget(
                os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_fpkm.csv")
            )
        return tasks

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
        return df

    @staticmethod
    def read_count_table() -> pd.DataFrame:
        """Read and verify a RNAseq FPRKM count table."""
        df_counts = pd.read_csv(ProbeConfig().encode_count_table, sep="\t")
        counts = df_counts[df_counts["gene_id"].str.contains("ENSG")].copy()
        counts["clean_id"] = counts["gene_id"].apply(lambda x: x.split(".")[0])
        return counts[constants.COUNTS_COLUMNS]

    def read_gtf_file(self):
        return pd.read_parquet(self.input()["gtf"].path)[constants.GTF_COLUMNS]

    def filter_gene_of_interest(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter FPKM table to exclude the gene of interest."""
        # Filter using provided EnsembleID directly
        if SequenceConfig().ensemble_id:
            self.logger.debug(
                f"Filtering directly using EmsembleID {SequenceConfig().ensemble_id}"
            )
            return df[df["clean_id"] != SequenceConfig().ensemble_id]

        # Filter using blast
        fname = os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_blast.txt")
        args_blast = [
            "blastn",
            "-db",
            util.get_genome_name(),
            "-query",
            self.input()["blastseq"].path,
            "-task",
            "megablast",
            "-outfmt",
            "6",
            "-num_threads",
            GeneralConfig().threads,
            ">",
            fname,
        ]
        subprocess.check_call(args_blast)
        df_blast = pd.read_csv(fname, sep="\t", names=constants.BLAST_COLUMNS)
        df_blast = df_blast[
            (df_blast["pident"] >= 98)
            & (df_blast["evalue"] <= 1e-8)
            & (df_blast["sseqid"] == df_blast.loc[0, "sseqid"])
        ]

        # Pad blast start/end results in case not everything was found/isoforms/etc.
        # Set to an arbitrary maximum value of 100 to not include other genes
        for buffer in range(0, 110, 10):
            gene_ids = []
            for _, row in df_blast.iterrows():
                gene_ids.extend(
                    df[
                        (df["seqname"] == row["sseqid"])
                        & (df["start"] >= row["sstart"] - buffer)
                        & (df["end"] <= row["send"] + buffer)
                    ]["clean_id"].values
                )
            if not gene_ids:
                self.logger.debug(f"No genes found using a buffer of {buffer}")
                continue

            most_frequent_gene_id = max(set(gene_ids), key=gene_ids.count)
            self.logger.debug(
                f"Filtering out most frequent gene {most_frequent_gene_id}"
            )
            df = df[df["clean_id"] != most_frequent_gene_id].reset_index(drop=True)
            break

        return df

    def get_maximum_fpkm(
        self, df_sam: pd.DataFrame, df_counts: pd.DataFrame, df_gtf: pd.DataFrame
    ) -> pd.DataFrame:
        """Filter candidates based on FPKM."""
        # Merge gtf and count table to get FPKM at each genomic location
        df_fpkm = pd.merge(
            df_counts,
            df_gtf,
            how="left",
            left_on="clean_id",
            right_on="gene_id",
        ).dropna()

        if self.is_blast_required:
            df_fpkm = self.filter_gene_of_interest(df_fpkm)

        # Match FPKM to the probe candidates' off-target genomic locations
        # Loop over chromosome/rname to keep loci separated
        dfs = []
        for sequence, df_sequence in tqdm(df_sam.groupby("rname")):
            df_fpkm_sequence = df_fpkm[df_fpkm["seqname"] == sequence]
            gene_start = df_fpkm_sequence["start"].astype(int).values
            gene_end = df_fpkm_sequence["end"].astype(int).values
            align_start = df_sequence["pos"].astype(int).values
            align_end = align_start + df_sequence["seq"].apply(len).values

            # Create n*m matrix of which alignments(n) match up with which annotations(m)
            i, j = np.where(
                (
                    (align_start[:, None] >= gene_start)
                    & (align_start[:, None] <= gene_end)
                )
                | (
                    (align_end[:, None] >= gene_start)
                    & (align_end[:, None] <= gene_end)
                )
            )

            # Join columns to merge alignment information with FPKM values
            df_sequence = pd.concat(
                [
                    df_sequence.reset_index(drop=True).iloc[i].reset_index(drop=True),
                    df_fpkm_sequence.reset_index(drop=True)
                    .iloc[j]
                    .reset_index(drop=True),
                ],
                axis=1,
            )
            dfs.append(df_sequence)
        dfs = pd.concat(dfs, ignore_index=True)

        # Get the highest FPKM per alignment/candidate
        df_max_fpkm = dfs.groupby("qname", as_index=False)["FPKM"].max()
        return df_max_fpkm

    def run(self):
        fname_fasta = self.input()["probes"].path
        fname_sam = os.path.splitext(fname_fasta)[0] + ".sam"
        self.align_probes(fname_fasta, fname_sam)
        df_sam = self.filter_unique_probes(fname_sam)

        if self.is_blast_required:
            df_gtf = self.read_gtf_file()
            df_counts = self.read_count_table()
            df_max_fpkm = self.get_maximum_fpkm(df_sam, df_counts, df_gtf)
            df_max_fpkm.to_csv(self.output()["table"].path, index=False)
            df_sam = df_max_fpkm[
                df_max_fpkm["FPKM"] <= ProbeConfig().max_off_target_fpkm
            ]

        # Filter candidates found in samfile
        sequences = list(Bio.SeqIO.parse(fname_fasta, "fasta"))
        candidates = [
            seq for seq in tqdm(sequences) if seq.id in df_sam["qname"].values
        ]
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, "fasta")
