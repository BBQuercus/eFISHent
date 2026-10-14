"""Filter probes based on alignment score and uniqueness."""

import logging
import os
import subprocess

import Bio.SeqIO
import luigi
import pandas as pd
import pysam

from . import constants, util
from .basic_filtering import BasicFiltering
from .config import GeneralConfig, ProbeConfig, SequenceConfig
from .indexing import BuildBowtie2Index, BuildBowtieIndex
from .off_target_annotation import AnnotationOffTargetFilterMixin
from .off_target_genome import GenomeOffTargetFilterMixin
from .rdna_filter import FilterRibosomalRNA


class AlignProbeCandidates(
    AnnotationOffTargetFilterMixin, GenomeOffTargetFilterMixin, luigi.Task
):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_endogenous = SequenceConfig().is_endogenous
        self.max_off_targets = ProbeConfig().max_off_targets
        self.no_alternative_loci = ProbeConfig().no_alternative_loci
        self.aligner = ProbeConfig().aligner
        self.encode_count_table = ProbeConfig().encode_count_table
        self.off_target_min_tm = ProbeConfig().off_target_min_tm
        self.mask_repeats = ProbeConfig().mask_repeats
        self.intergenic_off_targets = ProbeConfig().intergenic_off_targets
        self.filter_rrna = ProbeConfig().filter_rrna
        self.has_expression_filter = bool(
            self.encode_count_table
            and GeneralConfig().reference_annotation
            and self.is_endogenous
        )
        self.has_annotation = bool(GeneralConfig().reference_annotation)
        self.has_intergenic_filter = bool(
            self.intergenic_off_targets
            and self.has_annotation
            and self.is_endogenous
        )
        self.has_rrna_filter = bool(
            self.filter_rrna
            and self.has_annotation
        )
        self.has_pseudogene_filter = bool(
            self.has_annotation
            and self.is_endogenous
        )
        self.has_transcriptome = bool(GeneralConfig().reference_transcriptome)

    def requires(self):
        if self.aligner == "bowtie2":
            index_task = BuildBowtie2Index()
        else:
            index_task = BuildBowtieIndex()
        # FilterRibosomalRNA sits between BasicFiltering and alignment
        if ProbeConfig().filter_rdna_45s:
            probe_task = FilterRibosomalRNA()
        else:
            probe_task = BasicFiltering()
        tasks = {"probes": probe_task, "bowtie": index_task}
        if self.has_expression_filter or self.has_intergenic_filter or self.has_rrna_filter or self.has_pseudogene_filter:
            from .indexing import PrepareAnnotationFile
            tasks["annotation"] = PrepareAnnotationFile()
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

    def align_probes(self, threads: int):
        """Align probes to the reference genome.

        Takes a fasta file as input (which will be converted to fastq).
        Saves a sam file with alignments.
        """
        # Convert fasta to fastq - bowtie only allows uses fastq inputs
        fname_fastq = self.fname_fasta.rstrip("a") + "q"
        with open(self.fname_fasta) as fasta, open(fname_fastq, "w") as fastq:
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
        from .console import spinner

        with spinner("Aligning probes to genome..."):
            subprocess.check_call(
                args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )

    def align_probes_bowtie2(self, threads: int):
        """Align probes using Bowtie2 with OligoMiner/Tigerfish parameters.

        Parameters are from OligoMiner (Beliveau et al., 2018, MIT license)
        and validated by Tigerfish (Kuo et al., 2024). Tuned for short
        oligonucleotide probe alignment with mismatch tolerance.
        """
        # Bowtie2 misparses FASTA descriptions containing angle brackets
        # (e.g., BioPython's "<unknown description>"). Write a clean FASTA
        # with IDs only to avoid read name corruption.
        fname_clean = self.fname_fasta + ".bt2.fa"
        with open(fname_clean, "w") as out:
            for record in Bio.SeqIO.parse(self.fname_fasta, "fasta"):
                out.write(f">{record.id}\n{record.seq}\n")

        k_value = str(max(self.max_off_targets + 2, 100))

        # Adaptive seed length: use shorter seeds for short probes (<=22nt)
        # to avoid the blind spot where -L 20 makes Bowtie2 non-functional
        # (seed is essentially probe-length, so any mismatch causes a miss)
        min_probe_len = ProbeConfig().min_length
        if min_probe_len <= 22:
            seed_length = "10"
        else:
            seed_length = "20"

        if self.is_endogenous:
            # Local alignment for endogenous: detect partial off-target hits
            args = [
                "bowtie2",
                "-f",  # FASTA input
                "-x", self.fname_genome,
                "-U", fname_clean,
                "--threads", str(threads),
                "--local",  # Local alignment (soft-clip ends)
                "-D", "20",  # 20 seed extension attempts
                "-R", "3",  # 3 re-seeding rounds
                "-N", "1",  # 1 mismatch allowed in seed
                "-L", seed_length,  # Seed length (adaptive for short probes)
                "-i", "C,4",  # Seed interval: constant, every 4 bases
                "--score-min", "G,1,4",  # Min score: log-based threshold
                "-k", k_value,  # Report up to k alignments per probe
                "--no-unal",  # Suppress unaligned reads
                "-S", self.fname_sam,
            ]
        else:
            # End-to-end alignment for exogenous: only flag probes that
            # match the non-target genome across (nearly) their full length.
            # Local mode produces spurious partial hits that incorrectly
            # eliminate probes from a foreign sequence.
            args = [
                "bowtie2",
                "-f",  # FASTA input
                "-x", self.fname_genome,
                "-U", fname_clean,
                "--threads", str(threads),
                "--end-to-end",  # Require full-length alignment
                "-D", "20",  # 20 seed extension attempts
                "-R", "3",  # 3 re-seeding rounds
                "-N", "1",  # 1 mismatch allowed in seed
                "-L", seed_length,  # Seed length (adaptive for short probes)
                "-i", "C,4",  # Seed interval: constant, every 4 bases
                "--score-min", "L,-0.6,-0.6",  # Allow ~1 mismatch per 10 bases
                "-k", k_value,  # Report up to k alignments per probe
                "-S", self.fname_sam,
            ]
        self.logger.debug(f"Running bowtie2 with - {' '.join(args)}")
        from .console import spinner

        with spinner("Aligning probes to genome (Bowtie2)..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )
        os.remove(fname_clean)

    @staticmethod
    def _cigar_match_fraction(cigar: str) -> float:
        """Calculate fraction of query bases that are matched (M/=/X) in a CIGAR string.

        For bowtie2 local alignment, soft-clipped bases (S) indicate unaligned
        portions. Returns matched_bases / total_bases.
        """
        import re as _re
        ops = _re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if not ops:
            return 0.0
        total = 0
        matched = 0
        for length_str, op in ops:
            length = int(length_str)
            if op in ("M", "=", "X", "I", "S"):
                total += length
            if op in ("M", "=", "X"):
                matched += length
        return matched / total if total > 0 else 0.0

    @staticmethod
    def parse_raw_pysam(sam: str) -> pd.DataFrame:
        """Convert string based pysam output to a DataFrame."""
        # Parse tab and newline delimited pysam output
        data = [row.split("\t")[:10] for row in sam.split("\n")]
        columns = constants.SAMFILE_COLUMNS[:10]

        # If empty [['']]
        if data[0][0]:
            return pd.DataFrame(data, columns=columns).dropna()
        return pd.DataFrame(columns=columns)

    def filter_unique_probes(self) -> tuple:
        """Filter sam file with probes based on alignment score and uniqueness.

        Filtering explanations:
            * mapq 60: uniquely mapped read, regardless of number of mismatches / indels
            * flag 4: unmapped read
        For Bowtie2:
            * Endogenous: keep all mapped reads, filter by hit count
            * Exogenous: keep only unmapped (flag 4)

        Returns (df_filtered, df_pre_count) where df_pre_count is the DataFrame
        after mask/intergenic filtering but before count-based filtering (used by
        thermodynamic rescue to evaluate probes against the same filtered hits).
        """
        df = self.parse_raw_pysam(self._view_filtered_sam())
        df = self._drop_partial_alignments(df)
        self._warn_if_all_exogenous_unmapped(df)
        df = self._filter_pre_count_off_targets(df)

        # Snapshot before count-based filtering (for thermodynamic rescue)
        df_pre_count = df.copy()

        df = self._filter_by_off_target_count(df)
        return df, df_pre_count

    def _view_filtered_sam(self) -> str:
        """Return the SAM records relevant for filtering via samtools view."""
        if self.aligner != "bowtie2":
            flags = ["--min-MQ", "60"] if self.is_endogenous else ["--require-flags", "4"]
            self.logger.debug(f"Running samtools with - {' '.join(flags)}")
            return pysam.view(self.fname_sam, *flags)  # type: ignore
        if self.is_endogenous:
            # Keep all mapped reads, filter by hit count later
            return pysam.view(self.fname_sam, "--exclude-flags", "4")  # type: ignore
        return pysam.view(self.fname_sam, "--require-flags", "4")  # type: ignore

    def _drop_partial_alignments(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove bowtie2 local alignments that cover too little of the probe.

        The sensitive local parameters (-L 10) report matches of 10-11nt on
        20+nt probes — these won't stably hybridize. Keep only alignments
        where >=80% of the probe is matched.
        """
        if self.aligner != "bowtie2" or not self.is_endogenous or df.empty:
            return df

        pre_filter = len(df)
        min_match_frac = 0.80
        df = df[df["cigar"].apply(
            lambda c: self._cigar_match_fraction(c) >= min_match_frac
        )]
        n_removed = pre_filter - len(df)
        if n_removed > 0:
            self.logger.debug(
                f"Filtered {n_removed} partial alignments "
                f"(<80% match by CIGAR) from {pre_filter} total"
            )
        return df

    def _warn_if_all_exogenous_unmapped(self, df: pd.DataFrame) -> None:
        """Warn when all exogenous probes pass unmapped.

        This indicates alignment parameters may be too strict to find off-targets.
        """
        if self.is_endogenous or self.aligner != "bowtie2" or not hasattr(self, "fname_fasta"):
            return

        total_probes = len(list(Bio.SeqIO.parse(self.fname_fasta, "fasta")))
        if total_probes > 0 and len(df) == total_probes:
            self.logger.warning(
                f"All {total_probes} exogenous probes passed Bowtie2 unmapped. "
                "This may indicate alignment parameters are too strict to detect "
                "off-targets. The rDNA/satellite filter and transcriptome BLAST "
                "provide additional protection."
            )

    def _filter_pre_count_off_targets(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove hits that should not count towards a probe's off-targets."""
        if self.no_alternative_loci:
            df = df[~df["rname"].str.contains("_alt")]

        # Remove hits in repetitive/masked regions before counting
        if self.mask_repeats and self.is_endogenous:
            df = self.filter_masked_off_targets(df)

        # Remove hits in intergenic regions before counting
        if self.has_intergenic_filter:
            df = self.filter_intergenic_off_targets(df)

        # Remove hits overlapping pseudogene loci before counting
        # (pseudogenes produce little/no RNA — not real off-targets for RNA-FISH)
        if self.has_pseudogene_filter:
            df = self.filter_pseudogene_off_targets(df)

        return df

    def _filter_by_off_target_count(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove endogenous probes with more genomic hits than allowed."""
        if not self.is_endogenous:
            return df

        if self.has_transcriptome and self.max_off_targets == 0:
            # When a transcriptome is available and max_off_targets is
            # at the default (0), skip the genome count-based filter.
            # The downstream BLAST checks for expressed off-targets
            # directly, which is more accurate for RNA-FISH than raw
            # genome hit counts. Genome hits to pseudogenes and other
            # non-expressed regions inflate the count without reflecting
            # real off-target binding risk (e.g. ACTB has 100+ genomic
            # hits due to pseudogenes but only a handful of expressed
            # off-targets).
            self.logger.debug(
                "Skipping genome count-based off-target filter "
                "(transcriptome BLAST will handle off-target detection)"
            )
            return df

        group_sizes = df.groupby("qname").size()
        return df[
            df["qname"].isin(
                group_sizes[group_sizes <= self.max_off_targets + 1].index
            )
        ]

    def run(self):
        util.log_stage_start(self.logger, "AlignProbeCandidates")
        # Naming
        self.fname_fasta = self.input()["probes"].path
        self.fname_sam = os.path.splitext(self.fname_fasta)[0] + ".sam"
        self.fname_genome = util.get_genome_name()
        self.logger.info(
            f"Aligning with {self.aligner} against {os.path.basename(self.fname_genome)}"
        )

        # Alignment
        if self.aligner == "bowtie2":
            self.align_probes_bowtie2(threads=GeneralConfig().threads)
        else:
            self.align_probes(threads=GeneralConfig().threads)

        # Filtering
        df_sam, df_pre_count = self.filter_unique_probes()

        # Thermodynamic off-target filtering (rescue probes with unstable off-targets)
        if self.off_target_min_tm > 0 and self.is_endogenous:
            df_sam = self.filter_thermodynamic_off_targets(df_sam, df_pre_count)

        # rRNA off-target filtering
        if self.has_rrna_filter:
            df_sam = self.filter_rrna_off_targets(df_sam)

        # Expression-weighted filtering
        if self.has_expression_filter:
            df_sam = self.filter_expression_off_targets(df_sam)

        df_sam.to_csv(self.output()["table"].path, index=False)

        # Save output
        sequences = list(Bio.SeqIO.parse(self.fname_fasta, "fasta"))
        candidates = [seq for seq in sequences if seq.id in df_sam["qname"].unique()]
        self.logger.debug(df_sam["qname"].values)
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, format="fasta")
