"""
Create bowtie index for a given genome.
Align probes to the reference genome.
Filter probes based on alignment score and uniqueness.
"""

import logging
import os
import subprocess
import warnings

import Bio.SeqIO
import gtfparse
import luigi
import numpy as np
import pandas as pd
import pysam

from .config import GeneralConfig
from .config import ProbeConfig
from .config import SequenceConfig
from .basic_filtering import BasicFiltering
from .prepare_sequence import BuildBlastDatabase
from .prepare_sequence import PrepareSequence
from . import constants
from . import util


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
        ]
        self.logger.debug(f"Running bowtie with - {''.join(args_bowtie)}")
        subprocess.check_call(args_bowtie)


class PrepareAnnotationFile(luigi.Task):
    """Convert gtf annotation file to parquet.

    Doesn't change any data but greatly speeds up subsequent I/O (~10x increase).
    """

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(GeneralConfig().reference_annotation + ".parq")

    def prepare_gtf_file(self, fname_input: str, fname_output) -> None:
        """Save GTF file as parquet."""
        warnings.filterwarnings("ignore")
        df = gtfparse.read_gtf(fname_input)
        warnings.filterwarnings("default")
        self.logger.debug("Read gtf file with gtfparse")
        df.to_parquet(fname_output)

    def run(self):
        self.prepare_gtf_file(GeneralConfig().reference_annotation, self.output().path)


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

    def read_count_table(self, fname: str) -> pd.DataFrame:
        """Read and verify a RNAseq FPRKM count table."""
        df_counts = pd.read_csv(fname, sep="\t")
        counts = df_counts[df_counts["gene_id"].str.contains("ENSG")].copy()
        counts["clean_id"] = counts["gene_id"].apply(lambda x: x.split(".")[0])
        self.logger.debug(f"Read count table with {len(counts)} entries.")
        return counts[constants.COUNTS_COLUMNS]

    def read_gtf_file(self, fname: str):
        """Parse prepared parquet-formatted gtf file."""
        return pd.read_parquet(fname)[constants.GTF_COLUMNS]

    def align_probes(
        self,
        fname_input: str,
        fname_output: str,
        max_off_targets: int,
        is_endogenous: bool,
        threads: int,
    ):
        """Align probes to the reference genome.

        Takes a fasta file as input (which will be converted to fastq).
        Saves a sam file with alignments.
        """
        # Convert fasta to fastq - bowtie doesn't return read names if not in fastq...
        fname_fastq = fname_input.rstrip("a") + "q"
        with open(fname_input, "r") as fasta, open(fname_fastq, "w") as fastq:
            for record in Bio.SeqIO.parse(fasta, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                Bio.SeqIO.write(record, fastq, "fastq")
        self.logger.debug(f"Converted fasta to fastq - {fname_fastq}")

        alignments = str(max(max_off_targets + 1, 1))
        endo_exo_param = [] if is_endogenous else ["-m", alignments]
        args_bowtie = [
            "bowtie",
            util.get_genome_name(),
            fname_fastq,
            "--threads",
            str(threads),
            *endo_exo_param,
            "--sam",
            "-S",
            fname_output,
        ]
        self.logger.debug(f"Running bowtie with - {' '.join(args_bowtie)}")
        subprocess.check_call(
            args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )

    def filter_unique_probes(self, fname: str, is_endogenous: bool) -> pd.DataFrame:
        """Filter sam file with probes based on alignment score and uniqueness."""
        # 60 - uniquely mapped read, regardless of number of mismatches / indels
        # 4 – flag for unmapped read
        flags = ["--min-MQ", "60"] if is_endogenous else ["--require-flags", "4"]
        end = 14 if is_endogenous else 12
        self.logger.debug(f"Running samtools with - {' '.join(flags)}")
        filtered_sam = pysam.view(fname, *flags)

        # Parse tab and newline delimited pysam output
        df = pd.DataFrame(
            [row.split("\t") for row in filtered_sam.split("\n")],
            columns=constants.SAMFILE_COLUMNS[:end],
        )
        return df

    def filter_gene_of_interest(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter FPKM table to exclude the gene of interest."""
        # Filter using provided EnsemblID directly
        if SequenceConfig().ensembl_id:
            self.logger.debug(
                f"Filtering directly using EmsembleID {SequenceConfig().ensembl_id}"
            )
            return df[df["clean_id"] != SequenceConfig().ensembl_id]

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
        self.logger.debug(f"Running blast with - {''.join(args_blast)}")
        subprocess.check_call(args_blast)
        df_blast = pd.read_csv(fname, sep="\t", names=constants.BLAST_COLUMNS)

        # Arbirary identity values to ensure similar targets aren't falsely selected
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
        """Find the highest FPKM values across detected off-targets."""
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
        for sequence, df_sequence in df_sam.groupby("rname"):
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
        # Alignment
        fname_fasta = self.input()["probes"].path
        fname_sam = os.path.splitext(fname_fasta)[0] + ".sam"
        self.align_probes(
            fname_input=fname_fasta,
            fname_output=fname_sam,
            max_off_targets=ProbeConfig().max_off_targets,
            is_endogenous=SequenceConfig().is_endogenous,
            threads=GeneralConfig().threads,
        )
        df_sam = self.filter_unique_probes(
            fname=fname_sam, is_endogenous=SequenceConfig().is_endogenous
        )

        # TODO allow filtering if exogenous?
        # FPKM stuff with count table
        if self.is_blast_required:
            df_gtf = self.read_gtf_file(self.input()["gtf"].path)
            df_counts = self.read_count_table(ProbeConfig().encode_count_table)
            df_max_fpkm = self.get_maximum_fpkm(df_sam, df_counts, df_gtf)
            df_max_fpkm.to_csv(self.output()["table"].path, index=False)
            df_sam = df_max_fpkm[
                df_max_fpkm["FPKM"] <= ProbeConfig().max_off_target_fpkm
            ]

        # Filter candidates found in samfile
        sequences = list(Bio.SeqIO.parse(fname_fasta, "fasta"))
        candidates = [seq for seq in sequences if seq.id in df_sam["qname"].values]
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, "fasta")
