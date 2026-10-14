"""Clean up the output directory.

Remove the files that are not needed and prettify the kept output.
"""

import glob
import logging
import multiprocessing
import os
import sys

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import luigi
import pandas as pd

from . import util
from .alignment import AlignProbeCandidates
from .basic_filtering import (
    compute_duplex_dg,
    get_cpg_fraction,
    get_gc_content,
    get_melting_temp,
)
from .config import GeneralConfig, ProbeConfig, RunConfig, SequenceConfig
from .genome_verification import GenomeVerificationMixin
from .kmers import BuildJellyfishIndex, get_max_kmer_counts_batch
from .off_target_clustering import OffTargetClusteringMixin
from .off_target_report import OffTargetReportMixin
from .optimization import OptimizeProbeCoverage
from .probe_summary import ProbeSummaryMixin
from .quality import QualityScoringMixin
from .secondary_structure import get_free_energy


class CleanUpOutput(
    QualityScoringMixin,
    OffTargetReportMixin,
    OffTargetClusteringMixin,
    GenomeVerificationMixin,
    ProbeSummaryMixin,
    luigi.Task,
):
    """Clean up the output files and remove the intermediaries that are not needed."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_blast_required = (
            ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
        )

    def requires(self):
        tasks = {
            "optimize": OptimizeProbeCoverage(),
            "jellyfish": BuildJellyfishIndex(),
            "alignment": AlignProbeCandidates(),
        }
        return tasks

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("fasta", f"{util.get_gene_name()}.fasta"),
                ("table", f"{util.get_gene_name()}.csv"),
                ("config", f"{util.get_gene_name()}.txt"),
            ]
        }

    def prettify_table(
        self,
        sequences: list[Bio.SeqRecord.SeqRecord],
        basename: str,
        jellyfish_path: str,
        alignment_path: str,
    ) -> pd.DataFrame:
        """Create table with probe information."""
        cfg = ProbeConfig()

        # Create basic table with probe start/end positions and name
        df = util.create_data_table(sequences)
        sequences = sorted(sequences, key=lambda x: int(x.id.split("-")[-1]))

        # Add data columns
        df["GC"] = [round(get_gc_content(seq.seq), 2) for seq in sequences]
        df["TM"] = [
            round(
                get_melting_temp(
                    seq.seq, cfg.na_concentration, cfg.formamide_concentration
                ),
                2,
            )
            for seq in sequences
        ]
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            df["deltaG"] = pool.map(get_free_energy, sequences)
        df["kmers"] = get_max_kmer_counts_batch(sequences, jellyfish_path)
        df_counts = pd.read_csv(alignment_path)
        df_counts = df_counts.groupby("qname", as_index=False).size()
        df["count"] = pd.merge(
            df, df_counts, how="left", left_on="name", right_on="qname",
            validate="many_to_one",
        )["size"].fillna(0)
        df["count"] -= 1

        # Add CpG fraction and low-complexity score for quality scoring
        df["cpg_fraction"] = [
            round(get_cpg_fraction(seq.seq), 4) for seq in sequences
        ]
        df["low_complexity"] = [
            self._compute_low_complexity_score(seq.seq) for seq in sequences
        ]

        # Add target site accessibility scoring if enabled
        if cfg.accessibility_scoring:
            df["accessibility"] = self._compute_accessibility_scores(df)
        else:
            df["accessibility"] = 1.0

        # Add on-target duplex ΔG (binding stability)
        df["on_target_dg"] = [
            round(
                compute_duplex_dg(
                    seq.seq, cfg.na_concentration, cfg.formamide_concentration
                ),
                2,
            )
            for seq in sequences
        ]

        # Add transcriptome off-target details if available
        self._add_transcriptome_details(df)

        # Add expression risk annotation if available
        self._annotate_expression_risk(df)

        # Compute per-probe quality score (0-100)
        df["quality"] = self._compute_quality_scores(df)

        # Add recommendation column
        df["recommendation"] = df.apply(self._compute_recommendation, axis=1)

        # Create new/clean names
        df["name"] = [f"{basename}-{idx + 1}" for idx in df.index]
        return df

    def _compute_accessibility_scores(self, df: pd.DataFrame) -> pd.Series:
        """Compute target RNA accessibility for each probe's binding site.

        Uses local RNA folding (200nt window) around each probe binding site
        to determine what fraction of the target is unpaired (accessible).
        Returns a Series of accessibility scores in [0, 1].
        """
        from .accessibility import compute_accessibility_scores

        # Load the target sequence
        sequence_fasta = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_sequence.fasta"
        )
        if not os.path.isfile(sequence_fasta):
            cfg = SequenceConfig()
            if cfg.sequence_file and os.path.isfile(cfg.sequence_file):
                sequence_fasta = cfg.sequence_file
            else:
                self.logger.debug("No target sequence for accessibility scoring")
                return pd.Series(1.0, index=df.index)

        try:
            target_seq = str(next(Bio.SeqIO.parse(sequence_fasta, "fasta")).seq)
        except Exception:
            return pd.Series(1.0, index=df.index)

        return compute_accessibility_scores(df, target_seq)

    def prettify_sequences(self, df: pd.DataFrame) -> list[Bio.SeqRecord.SeqRecord]:
        """Clean up sequence id's and descriptions."""
        return [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(seq), id=name, name=name, description=""
            )
            for seq, name in zip(df["sequence"], df["name"])
        ]

    def prettify_section(self, section: str) -> str:
        """Prettify a single configuration section."""
        pretty = [f"{section.replace('Config', '')} configuration:"]
        for key, value in self.config[section].items():
            if value and value != '""':
                pretty.append(f"  {key.replace('_', '-')}: {value}")
        return "\n".join(pretty)

    def prettify_configuration(self) -> str:
        """Create a pretty file with all set/default configuration."""
        self.config = luigi.configuration.get_config()
        pretty = [self.prettify_section(section) for section in self.config.sections()]
        command = " ".join(sys.argv[1:])
        pretty.append(f"Command:\n  efishent {command}")
        return "\n\n".join(pretty)

    def remove_intermediates(self) -> None:
        """Remove intermediary files no longer used."""
        intermediary_files = glob.glob(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_*")
        )
        self.logger.debug(f"Found file intermediates - {intermediary_files}")
        try:
            intermediary_files.remove(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_entrez.fasta"
                )
            )
        except ValueError:
            self.logger.debug("Didn't find entrez file...")
        for filename in intermediary_files:
            self.logger.debug(f"Removing {filename}")
            os.remove(filename)

    def run(self):
        util.log_stage_start(self.logger, "CleanUpOutput")
        sequences = list(
            Bio.SeqIO.parse(self.input()["optimize"]["probes"].path, "fasta")  # type: ignore
        )
        df = self.prettify_table(
            sequences,
            basename=util.get_gene_name(hashed=False),
            jellyfish_path=self.input()["jellyfish"].path,  # type: ignore
            alignment_path=self.input()["alignment"]["table"].path,  # type: ignore
        )

        # Apply cumulative off-target cap (C3)
        df = self._apply_off_target_cap(df)

        # Detect off-target clustering and emit warnings
        clustering_warnings = self._detect_off_target_clustering(df)
        if clustering_warnings:
            from .console import print_warning
            for warning in clustering_warnings:
                self.logger.warning(warning)
                print_warning(warning)

        sequences = self.prettify_sequences(df)
        config = self.prettify_configuration()

        # Write initial output for BLAST verification
        Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")

        # Run BLAST verification and remove flagged probes
        self.logger.info("Verifying probes against genome (BLAST)")
        self._verification = self._run_blast_verification(
            self.output()["fasta"].path
        )
        if self._verification and self._verification["flagged"]:
            flagged_names = set(self._verification["flagged"].keys())
            pre_count = len(df)
            df = df[~df["name"].isin(flagged_names)].reset_index(drop=True)
            self.logger.debug(
                f"BLAST verification removed {pre_count - len(df)} probes with "
                f"unexpected genomic hits: {sorted(flagged_names)[:5]}"
                f"{'...' if len(flagged_names) > 5 else ''}"
            )
            # Rewrite output with cleaned probes
            sequences = self.prettify_sequences(df)
            Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")
            # Update verification counts
            self._verification["total"] = len(df)
            self._verification["clean"] = len(df)
            self._verification["flagged"] = {}

        util.log_and_check_candidates(self.logger, "CleanUpOutput", len(sequences))

        # Minimum probe count warning
        min_probes = 25
        if 0 < len(sequences) < min_probes:
            from .console import is_silent, print_warning
            avg_quality = df["quality"].mean() if "quality" in df.columns else 0
            gc_range = f"{df['GC'].min():.0f}-{df['GC'].max():.0f}%" if "GC" in df.columns else "N/A"
            msg = (
                f"Low probe count: only {len(sequences)} probes in final set "
                f"(recommended >=25 for smFISH). "
                f"Mean quality: {avg_quality:.0f}/100, GC range: {gc_range}. "
                f"Note: probe quality matters more than quantity — "
                f"high-quality 12-probe sets can outperform low-quality 70-probe sets. "
                f"To get more probes: try relaxing --min-tm/--max-gc, "
                f"reducing --spacing, or using --target-regions both."
            )
            self.logger.warning(msg)
            if not is_silent():
                print_warning(msg)

        # Store summary stats and probe data for completion message
        self._summary = self._compute_summary(df)
        self._probe_df = df

        df.to_csv(self.output()["table"].path, index=False)
        with open(self.output()["config"].path, "w") as f:
            f.write(config)
        self.logger.debug(f'Saving all files using hash "{util.get_gene_name()}"')

        # Files to be deleted
        if not RunConfig().save_intermediates:
            self.remove_intermediates()
