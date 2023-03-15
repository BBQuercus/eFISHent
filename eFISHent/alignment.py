"""Filter probes based on alignment score and uniqueness."""

# from typing import Any
import logging
import os
import subprocess

# import numpy as np
import Bio.SeqIO
import luigi
import pandas as pd
import pysam

# from .indexing import BuildBlastDatabase
# from .indexing import PrepareAnnotationFile
# from .prepare_sequence import PrepareSequence
from . import constants
from . import util
from .basic_filtering import BasicFiltering
from .config import GeneralConfig
from .config import ProbeConfig
from .config import SequenceConfig
from .indexing import BuildBowtieIndex


class AlignProbeCandidates(luigi.Task):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.is_blast_required = (
        #     ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
        # )
        self.is_endogenous = SequenceConfig().is_endogenous
        self.max_off_targets = ProbeConfig().max_off_targets

    def requires(self):
        tasks = {"probes": BasicFiltering(), "bowtie": BuildBowtieIndex()}
        # if self.is_blast_required:
        #     tasks["blastseq"] = PrepareSequence()
        #     tasks["gtf"] = PrepareAnnotationFile()
        #     tasks["blast"] = BuildBlastDatabase()
        return tasks

    def output(self):
        tasks = {
            "fasta": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_aligned.fasta"
                )
            ),
            "table": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_counts.csv"
                )
            ),
        }
        return tasks

    @property
    def count_table(self) -> pd.DataFrame:
        """Read and verify a RNAseq normalized count table."""
        extension = os.path.splitext(self.fname_count)[1].lower()
        if extension not in (".tsv", ".csv", ".txt"):
            raise ValueError(
                "The count table must be provided as TSV, CSV, or TXT file. "
                f"Only found - {self.fname_count}"
            )
        sep = "\t" if extension in (".tsv", ".txt") else ","
        df = pd.read_csv(self.fname_count, sep=sep, index_col=0).reset_index().dropna()
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

    def align_probes(self, threads: int):
        """Align probes to the reference genome.

        Takes a fasta file as input (which will be converted to fastq).
        Saves a sam file with alignments.
        """
        # Convert fasta to fastq - bowtie only allows uses fastq inputs
        fname_fastq = self.fname_fasta.rstrip("a") + "q"
        with open(self.fname_fasta, "r") as fasta, open(fname_fastq, "w") as fastq:
            for record in Bio.SeqIO.parse(fasta, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                Bio.SeqIO.write(record, fastq, format="fastq")
        self.logger.debug(f"Converted fasta to fastq - {fname_fastq}")

        # Actual alignment with -m to only return alignments with >m hits
        alignments = str(max(self.max_off_targets + 1, 1))
        endo_exo_param = [] if self.is_endogenous else ["-m", alignments]
        args_bowtie = [
            "bowtie",
            fname_fastq,
            "-x",
            self.fname_genome,
            "--threads",
            str(threads),
            *endo_exo_param,
            "--sam",
            self.fname_sam,
            "--all",
        ]
        self.logger.debug(f"Running bowtie with - {' '.join(args_bowtie)}")
        subprocess.check_call(
            args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )

    @staticmethod
    def parse_raw_pysam(sam: str) -> pd.DataFrame:
        """Convert string based pysam output to a DataFrame."""
        # Parse tab and newline delimited pysam output
        data = [row.split("\t")[:10] for row in sam.split("\n")]
        # end = 14 if is_endogenous else 12
        columns = constants.SAMFILE_COLUMNS[:10]

        # If empty [['']]
        if data[0][0]:
            return pd.DataFrame(data, columns=columns).dropna()
        return pd.DataFrame(columns=columns)

    def filter_unique_probes(self) -> pd.DataFrame:
        """Filter sam file with probes based on alignment score and uniqueness."""
        # 60 - uniquely mapped read, regardless of number of mismatches / indels
        # 4 – flag for unmapped read
        flags = ["--min-MQ", "60"] if self.is_endogenous else ["--require-flags", "4"]
        self.logger.debug(f"Running samtools with - {' '.join(flags)}")
        filtered_sam = pysam.view(self.fname_sam, *flags)  # type: ignore
        df = self.parse_raw_pysam(filtered_sam)

        # Filter to remove off-target rates
        if self.is_endogenous:
            group_sizes = df.groupby("qname").size()
            df = df[
                df["qname"].isin(
                    group_sizes[group_sizes <= self.max_off_targets + 1].index
                )
            ]
        return df

    # def exclude_gene_of_interest(
    #     self, df: pd.DataFrame, ensembl_id: str, fname_full_gene: str, threads: int
    # ) -> pd.DataFrame:
    #     """Filter count table to exclude the gene of interest."""
    #     # Filter using provided EnsemblID directly
    #     if ensembl_id:
    #         self.logger.debug(f"Filtering directly using EmsembleID {ensembl_id}")
    #         return df[df["gene_id"] != ensembl_id]

    #     # Filter using blast
    #     args_blast = [
    #         "blastn",
    #         "-db",
    #         self.fname_genome,
    #         "-query",
    #         fname_full_gene,
    #         "-task",
    #         "megablast",
    #         "-outfmt",
    #         "6",
    #         "-num_threads",
    #         str(threads),
    #     ]
    #     self.logger.debug(f"Running blast with - {' '.join(args_blast)}")
    #     process = subprocess.run(args_blast, capture_output=True, text=True)
    #     data = [row.split("\t") for row in process.stdout.split("\n") if row]
    #     df_blast = pd.DataFrame(data, columns=constants.BLAST_COLUMNS).astype(
    #         {name: float for name in ["pident", "sstart", "send", "evalue"]}
    #     )

    #     # Arbirary identity values to ensure similar targets aren't falsely selected
    #     df_blast = df_blast[
    #         (df_blast["pident"] >= 98)
    #         & (df_blast["evalue"] <= 1e-5)
    #         & (df_blast["sseqid"] == df_blast.loc[0, "sseqid"])
    #     ]

    #     # Pad blast start/end results in case not everything was found/isoforms/etc.
    #     # Set to an arbitrary maximum value of 100 to not include other genes
    #     # Account for forward and reverse direction
    #     for buffer in range(0, 110, 10):
    #         gene_ids = []
    #         for _, row in df_blast.iterrows():
    #             gene_ids.extend(
    #                 df[
    #                     (
    #                         (df["seqname"] == row["sseqid"])
    #                         | (df["seqname"] == f"chr{row['sseqid']}")
    #                     )
    #                     & ((df["start"] >= row["sstart"]) & (df["end"] <= row["send"]))
    #                     | ((df["start"] >= row["send"]) & (df["end"] <= row["sstart"]))
    #                 ]["gene_id"].values
    #             )
    #         if not gene_ids:
    #             self.logger.debug(f"No genes found using a buffer of {buffer}")
    #             continue

    #         most_frequent_gene_id = max(set(gene_ids), key=gene_ids.count)
    #         self.logger.debug(
    #             f"Filtering out most frequent gene {most_frequent_gene_id}"
    #         )
    #         df = df[df["gene_id"] != most_frequent_gene_id].reset_index(drop=True)
    #         break

    #     return df

    # def join_alignment_with_annotation(
    #     self, df_sam: pd.DataFrame, df_norm: pd.DataFrame
    # ) -> pd.DataFrame:
    #     """Merge probe alignment with normalized count table based on start/end positions."""
    #     dfs = []

    #     for sequence, df_sequence in df_sam.groupby("rname"):
    #         df_count_sequence = df_norm[
    #             (df_norm["seqname"] == sequence)
    #             | (df_norm["seqname"] == sequence.lstrip("chr"))
    #         ]
    #         gene_start = df_count_sequence["start"].astype(int).values
    #         gene_end = df_count_sequence["end"].astype(int).values
    #         align_start = df_sequence["pos"].astype(int).values
    #         align_length = df_sequence["seq"].apply(len).values
    #         align_strand = (
    #             (df_sequence["flag"] == constants.SAM_FLAG_REVERSE)
    #             .replace({True: -1, False: 1})
    #             .values
    #         )
    #         align_end = align_start + align_strand * align_length

    #         # Create n*m matrix of which alignments(n) match up with which annotations(m)
    #         i, j = np.where(
    #             (
    #                 (align_start[:, None] >= gene_start)
    #                 & (align_start[:, None] <= gene_end)
    #             )
    #             | (
    #                 (align_end[:, None] >= gene_start)
    #                 & (align_end[:, None] <= gene_end)
    #             )
    #         )

    #         # Join columns to merge alignment information with count values
    #         df_sequence = pd.concat(
    #             [
    #                 df_sequence.reset_index(drop=True).iloc[i].reset_index(drop=True),
    #                 df_count_sequence.reset_index(drop=True)
    #                 .iloc[j]
    #                 .reset_index(drop=True),
    #             ],
    #             axis=1,
    #         )
    #         dfs.append(df_sequence)
    #     dfs = pd.concat(dfs, ignore_index=True)
    #     return dfs

    # def get_most_expressed_genes(self, df: pd.DataFrame, percentage: float) -> Any:
    #     """Find the most highly expressed genes based on normalized count."""
    #     quantile = df["count"].quantile(percentage / 100)
    #     most_expressed = df[df["count"] > quantile].reset_index(drop=True)
    #     return most_expressed["gene_id"].values.tolist()

    # # TODO allow filtering if exogenous?
    # def filter_using_count(self, df_sam: pd.DataFrame) -> pd.DataFrame:
    #     """Filter candidates binding off targets above a expression percentage threshold."""
    #     df = self.join_alignment_with_annotation(df_sam=df_sam, df_norm=self.norm_table)
    #     df.to_csv(self.output()["table"].path, index=False)
    #     df_off_targets = self.exclude_gene_of_interest(
    #         df,
    #         ensembl_id=SequenceConfig().ensembl_id,
    #         fname_full_gene=self.input()["blastseq"].path,
    #         threads=GeneralConfig().threads,
    #     )
    #     most_expressed = self.get_most_expressed_genes(
    #         df=df_off_targets, percentage=ProbeConfig().max_expression_percentage
    #     )
    #     df_sam = df[~df["gene_id"].isin(most_expressed)].reset_index(drop=True)
    #     return df_sam

    def run(self):
        # Naming
        self.fname_fasta = self.input()["probes"].path
        self.fname_sam = os.path.splitext(self.fname_fasta)[0] + ".sam"
        self.fname_genome = util.get_genome_name()

        # Alignment and filtering
        self.align_probes(threads=GeneralConfig().threads)
        df_sam = self.filter_unique_probes()
        df_sam.to_csv(self.output()["table"].path, index=False)

        # if self.is_blast_required:
        #     self.fname_gtf = self.input()["gtf"].path
        #     self.fname_count = ProbeConfig().encode_count_table
        #     df_sam = self.filter_using_count(df_sam)

        # Save output
        sequences = list(Bio.SeqIO.parse(self.fname_fasta, "fasta"))
        candidates = [seq for seq in sequences if seq.id in df_sam["qname"].unique()]
        self.logger.debug(df_sam["qname"].values)
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, format="fasta")
