"""
Create bowtie index for a given genome.
Align probes to the reference genome.
Filter probes based on alignment score and uniqueness.
"""

import logging
import os
import subprocess
import tempfile
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
        subprocess.check_call(
            args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )
        self.logger.info("Finished building bowtie index.")


class PrepareAnnotationFile(luigi.Task):
    """Convert gtf annotation file to parquet.

    Doesn't change any data but greatly speeds up subsequent I/O (~10x increase).
    """

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(GeneralConfig().reference_annotation + ".parq")

    def prepare_gtf_file(self, fname_input: str, fname_output: str) -> None:
        """Save GTF file as parquet."""
        warnings.filterwarnings("ignore")
        df = gtfparse.read_gtf(fname_input)
        warnings.filterwarnings("default")
        self.logger.debug("Read gtf file with gtfparse")
        df.to_parquet(fname_output)

    def run(self):
        self.prepare_gtf_file(GeneralConfig().reference_annotation, self.output().path)
        self.logger.info("Finished parsing GTF annotation file.")


class AlignProbeCandidates(luigi.Task):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_blast_required = (
            ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
        )

    def requires(self):
        tasks = {"probes": BasicFiltering(), "bowtie": BuildBowtieIndex()}
        if self.is_blast_required:
            tasks["blastseq"] = PrepareSequence()
            tasks["gtf"] = PrepareAnnotationFile()
            tasks["blast"] = BuildBlastDatabase()
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
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_counts.csv"
                )
            )
        return tasks

    @property
    def count_table(self) -> pd.DataFrame:
        """Read and verify a RNAseq normalized count table."""
        if not os.path.splitext(self.fname_count)[1] == ".tsv":
            raise ValueError(
                "The count table must be provided as TSV file. "
                f"Only found - {self.fname_count}"
            )
        df = pd.read_csv(self.fname_count, sep="\t").reset_index().dropna()
        if not len(df.columns) >= 2:
            raise ValueError(
                "The count table must contain at least 2 columns."
                f"Only found - {df.columns}"
            )
        self.logger.debug(f"Mapping original the first two columns of - {df.columns}.")
        self.logger.debug(f"Read count table with {len(df)} entries.")
        self.logger.debug(f"Found dtypes - {df.dtypes}")
        self.logger.debug(f"First column - {df.iloc[0]}")

        df = df.iloc[:, :2]
        df.columns = ["gene_id", "count"]
        try:
            df = df.astype({"gene_id": str, "count": float})
        except ValueError:
            raise ValueError(
                "Could not convert count table columns to the right type. "
                "Please ensure the first two colums are of type string and number/float/int, respectively. "
            )

        # Remove potential version numbers
        df["gene_id"] = df["gene_id"].apply(lambda x: x.split(".")[0])
        return df

    @property
    def gtf_table(self) -> pd.DataFrame:
        """Parse prepared parquet-formatted gtf file."""
        return pd.read_parquet(self.fname_gtf)[constants.GTF_COLUMNS]

    @property
    def norm_table(self) -> pd.DataFrame:
        """Merge gtf and count table."""
        df = pd.merge(
            self.count_table, self.gtf_table, how="left", on="gene_id"
        ).dropna()
        return df

    def align_probes(self, max_off_targets: int, is_endogenous: bool, threads: int):
        """Align probes to the reference genome.

        Takes a fasta file as input (which will be converted to fastq).
        Saves a sam file with alignments.
        """
        # Convert fasta to fastq - bowtie doesn't return read names if not in fastq...
        fname_fastq = self.fname_fasta.rstrip("a") + "q"
        with open(self.fname_fasta, "r") as fasta, open(fname_fastq, "w") as fastq:
            for record in Bio.SeqIO.parse(fasta, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                Bio.SeqIO.write(record, fastq, "fastq")
        self.logger.debug(f"Converted fasta to fastq - {fname_fastq}")

        # Actual alignment with -m to only return alignments with >m hits
        alignments = str(max(max_off_targets + 1, 1))
        endo_exo_param = [] if is_endogenous else ["-m", alignments]
        args_bowtie = [
            "bowtie",
            self.fname_genome,
            fname_fastq,
            "--threads",
            str(threads),
            *endo_exo_param,
            "--sam",
            "-S",
            self.fname_sam,
        ]
        self.logger.debug(f"Running bowtie with - {' '.join(args_bowtie)}")
        subprocess.check_call(
            args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )

    def filter_unique_probes(self, is_endogenous: bool) -> pd.DataFrame:
        """Filter sam file with probes based on alignment score and uniqueness."""
        # 60 - uniquely mapped read, regardless of number of mismatches / indels
        # 4 – flag for unmapped read
        flags = ["--min-MQ", "60"] if is_endogenous else ["--require-flags", "4"]
        self.logger.debug(f"Running samtools with - {' '.join(flags)}")
        filtered_sam = pysam.view(self.fname_sam, *flags)

        # Parse tab and newline delimited pysam output
        data = [row.split("\t") for row in filtered_sam.split("\n")]
        end = 14 if is_endogenous else 12
        columns = constants.SAMFILE_COLUMNS[:end]

        # Padding in case some samfiles are larger, first cols will stay the same
        columns.extend(["" for _ in range(len(columns) - len(data[0]))])

        # If empty [['']]
        if data[0][0]:
            return pd.DataFrame(data, columns=columns).dropna()
        return pd.DataFrame(columns=columns)

    def exclude_gene_of_interest(
        self, df: pd.DataFrame, ensembl_id: str, fname_full_gene: str, threads: int
    ) -> pd.DataFrame:
        """Filter count table to exclude the gene of interest."""
        # Filter using provided EnsemblID directly
        if ensembl_id:
            self.logger.debug(f"Filtering directly using EmsembleID {ensembl_id}")
            return df[df["gene_id"] != ensembl_id]

        # Filter using blast
        args_blast = [
            "blastn",
            "-db",
            self.fname_genome,
            "-query",
            fname_full_gene,
            "-task",
            "megablast",
            "-outfmt",
            "6",
            "-num_threads",
            str(threads),
        ]
        self.logger.debug(f"Running blast with - {' '.join(args_blast)}")
        process = subprocess.run(args_blast, capture_output=True, text=True)
        data = [row.split("\t") for row in process.stdout.split("\n") if row]
        df_blast = pd.DataFrame(data, columns=constants.BLAST_COLUMNS).astype(
            {name: float for name in ["pident", "sstart", "send", "evalue"]}
        )

        # Arbirary identity values to ensure similar targets aren't falsely selected
        df_blast = df_blast[
            (df_blast["pident"] >= 98)
            & (df_blast["evalue"] <= 1e-5)
            & (df_blast["sseqid"] == df_blast.loc[0, "sseqid"])
        ]

        # Pad blast start/end results in case not everything was found/isoforms/etc.
        # Set to an arbitrary maximum value of 100 to not include other genes
        # Account for forward and reverse direction
        for buffer in range(0, 110, 10):
            gene_ids = []
            for _, row in df_blast.iterrows():
                gene_ids.extend(
                    df[
                        (
                            (df["seqname"] == row["sseqid"])
                            | (df["seqname"] == f"chr{row['sseqid']}")
                        )
                        & ((df["start"] >= row["sstart"]) & (df["end"] <= row["send"]))
                        | ((df["start"] >= row["send"]) & (df["end"] <= row["sstart"]))
                    ]["gene_id"].values
                )
            if not gene_ids:
                self.logger.debug(f"No genes found using a buffer of {buffer}")
                continue

            most_frequent_gene_id = max(set(gene_ids), key=gene_ids.count)
            self.logger.debug(
                f"Filtering out most frequent gene {most_frequent_gene_id}"
            )
            df = df[df["gene_id"] != most_frequent_gene_id].reset_index(drop=True)
            break

        return df

    def join_alignment_with_annotation(
        self, df_sam: pd.DataFrame, df_norm: pd.DataFrame
    ) -> pd.DataFrame:
        """Merge probe alignment with normalized count table based on start/end positions."""
        dfs = []

        for sequence, df_sequence in df_sam.groupby("rname"):
            df_count_sequence = df_norm[
                (df_norm["seqname"] == sequence)
                | (df_norm["seqname"] == sequence.lstrip("chr"))
            ]
            gene_start = df_count_sequence["start"].astype(int).values
            gene_end = df_count_sequence["end"].astype(int).values
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

            # Join columns to merge alignment information with count values
            df_sequence = pd.concat(
                [
                    df_sequence.reset_index(drop=True).iloc[i].reset_index(drop=True),
                    df_count_sequence.reset_index(drop=True)
                    .iloc[j]
                    .reset_index(drop=True),
                ],
                axis=1,
            )
            dfs.append(df_sequence)
        dfs = pd.concat(dfs, ignore_index=True)
        return dfs

    def get_most_expressed_genes(self, df: pd.DataFrame, percentage: float) -> list:
        """Find the most highly expressed genes based on normalized count."""
        quantile = df["count"].quantile(percentage / 100)
        most_expressed = df[df["count"] > quantile].reset_index(drop=True)
        return most_expressed["gene_id"].values.tolist()

    # TODO allow filtering if exogenous?
    def filter_using_count(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Filter candidates binding off targets with too high FPKMs."""
        df = self.exclude_gene_of_interest(
            self.norm_table,
            ensembl_id=SequenceConfig().ensembl_id,
            fname_full_gene=self.input()["blastseq"].path,
            threads=GeneralConfig().threads,
        )
        df = self.join_alignment_with_annotation(df_sam=df_sam, df_norm=df)
        most_expressed = self.get_most_expressed_genes(
            df=self.norm_table, percentage=ProbeConfig().max_expression_percentage
        )
        df.to_csv(self.output()["table"].path, index=False)
        df_sam = df[~df["gene_id"].isin(most_expressed)].reset_index(drop=True)
        return df_sam

    def run(self):
        self.logger.warning(f"is required - {self.is_blast_required}")
        self.logger.warning(ProbeConfig())

        # Naming
        self.fname_fasta = self.input()["probes"].path
        self.fname_sam = os.path.splitext(self.fname_fasta)[0] + ".sam"
        self.fname_genome = util.get_genome_name()
        self.fname_gene = util.get_gene_name()
        self.fname_gtf = self.input()["gtf"].path
        self.fname_count = ProbeConfig().encode_count_table

        # Alignment and filtering
        self.align_probes(
            max_off_targets=ProbeConfig().max_off_targets,
            is_endogenous=SequenceConfig().is_endogenous,
            threads=GeneralConfig().threads,
        )
        df_sam = self.filter_unique_probes(is_endogenous=SequenceConfig().is_endogenous)
        if self.is_blast_required:
            df_sam = self.filter_using_count(df_sam)

        # Save output
        sequences = list(Bio.SeqIO.parse(self.fname_fasta, "fasta"))
        candidates = [seq for seq in sequences if seq.id in df_sam["qname"].values]
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, "fasta")
