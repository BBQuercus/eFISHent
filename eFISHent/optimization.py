"""Find the optimal probeset with the highest coverage (fast or optimal)."""

import functools
import logging
import multiprocessing
import os

import Bio.Seq
import Bio.SeqIO
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import util
from .accessibility import compute_accessibility_scores, fold_sequence
from .config import GeneralConfig, ProbeConfig, RunConfig, SequenceConfig
from .secondary_structure import SecondaryStructureFiltering
from .selection import (
    cleanliness,
    compute_pre_optimization_quality,
    fill_coverage_gaps,
    greedy_model,
    is_overlapping,
    optimal_model,
    selection_weights,
)
from .tm_refinement import TmUniformityMixin

__all__ = [
    "OptimizeProbeCoverage",
    "cleanliness",
    "compute_accessibility_scores",
    "compute_pre_optimization_quality",
    "fill_coverage_gaps",
    "fold_sequence",
    "greedy_model",
    "is_binding",
    "is_overlapping",
    "optimal_model",
    "selection_weights",
    "visualize_assignment",
]


class OptimizeProbeCoverage(TmUniformityMixin, luigi.Task):
    """Find the optimal probeset with the highest coverage (fast or optimal)."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        if GeneralConfig().reference_transcriptome:
            from .transcriptome_filter import TranscriptomeFiltering
            return TranscriptomeFiltering()
        return SecondaryStructureFiltering()

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("probes", f"{util.get_gene_name()}_optimal.fasta"),
                ("coverage", f"{util.get_gene_name()}_optimization.png"),
            ]
        }

    def run_optimal(self, threads: int, time_limit: int):
        """Run the optimal model by parallelizing non-overlapping segments."""
        index_block = 0
        prev_probe = self.df.iloc[0]
        blocks = []

        for _, probe in self.df.iterrows():
            if not is_overlapping(
                (prev_probe["start"], prev_probe["end"]),
                (probe["start"], probe["end"]),
            ):
                index_block += 1
            prev_probe = probe
            blocks.append(index_block)
        self.df["block"] = blocks

        self.time_limit = time_limit
        with multiprocessing.Pool(threads) as pool:
            _assigned = pool.map(self.run_optimal_block, self.df["block"].unique())

        # Flatten the list of lists
        assigned = functools.reduce(lambda x, y: x + y, _assigned)
        return assigned

    def run_optimal_block(self, block: int) -> list[str]:
        """Run a single non-overlapping block."""
        df = self.df[self.df["block"] == block]
        return optimal_model(df, self.time_limit)

    def run_greedy(self):
        """Run sequential greedy model."""
        return greedy_model(self.df)

    def screen_genome_hits(
        self, sequences: list[Bio.SeqRecord.SeqRecord], fasta_path: str
    ) -> list[Bio.SeqRecord.SeqRecord]:
        """Drop exogenous candidates that the final genome BLAST check rejects.

        CleanUpOutput removes final probes with any near-full-length genomic
        hit. Running that same check on all candidates here lets the optimizer
        fill those positions with clean probes instead of leaving gaps.
        """
        verification = self._genome_verification(fasta_path)
        if not verification or not verification["flagged"]:
            return sequences
        flagged = set(verification["flagged"])
        if len(flagged) == len(sequences):
            self.logger.warning(
                "Every candidate has a genomic BLAST hit; keeping them for optimization."
            )
            return sequences
        self.logger.info(
            f"Genome BLAST screen removed {len(flagged)}/{len(sequences)} "
            "candidates with genomic hits"
        )
        return [seq for seq in sequences if seq.id not in flagged]

    @staticmethod
    def _genome_verification(fasta_path: str) -> dict | None:
        """CleanUpOutput's genome BLAST check on a candidate FASTA."""
        from .cleanup import CleanUpOutput

        return CleanUpOutput()._run_blast_verification(fasta_path)

    def filter_binding_probes(
        self, assigned: list[str], match_percentage: float
    ) -> list[str]:
        """Remove probes of too high rev/comp sequence similarity.

        Uses k-mer overlap (O(n²) comparisons but each is O(L) instead of
        O(L²) alignment) and upper-triangle-only to halve comparisons.
        """
        assigned_rows = self.df[self.df["name"].isin(assigned)]
        assigned_sequences = assigned_rows["sequence"].values
        n = len(assigned_sequences)

        # Upper triangle only — binding is symmetric
        trouble_indices = set()
        for i in range(n):
            for j in range(i + 1, n):
                if is_binding(assigned_sequences[i], assigned_sequences[j], match_percentage):
                    trouble_indices.add(j)

        self.logger.debug(f"Found {len(trouble_indices)} trouble makers to be removed.")
        assigned = [
            name for idx, name in enumerate(assigned) if idx not in trouble_indices
        ]
        return assigned

    def run(self):
        util.log_stage_start(self.logger, "OptimizeProbeCoverage")
        inp = self.input()
        fasta_path = inp["fasta"].path if isinstance(inp, dict) else inp.path
        sequences = list(Bio.SeqIO.parse(fasta_path, "fasta"))
        if not SequenceConfig().is_endogenous:
            sequences = self.screen_genome_hits(sequences, fasta_path)
        self.df = util.create_data_table(sequences)
        self.df["end"] += ProbeConfig().spacing

        # Off-target burden from TranscriptomeFiltering ranks otherwise-equal probes
        burden_target = inp.get("burden") if isinstance(inp, dict) else None
        if burden_target is not None and burden_target.exists():
            burden = pd.read_csv(burden_target.path).set_index("name")["burden"]
            self.df["_burden"] = self.df["name"].map(burden).fillna(0.0)

        # Compute accessibility scores for optimization weighting (not filtering)
        if ProbeConfig().accessibility_scoring:
            target_seq = self._load_target_sequence()
            if target_seq:
                self.logger.debug("Computing accessibility scores for optimization...")
                self.df["_accessibility"] = compute_accessibility_scores(
                    self.df, target_seq
                )
                n_accessible = (self.df["_accessibility"] > 0.5).sum()
                self.logger.debug(
                    f"Accessibility: {n_accessible}/{len(self.df)} probes "
                    f"target >50% accessible sites"
                )

        if RunConfig().optimization_method == "greedy":
            assigned = self.run_greedy()
        elif RunConfig().optimization_method == "optimal":
            assigned = self.run_optimal(
                threads=GeneralConfig().threads,
                time_limit=RunConfig().optimization_time_limit,
            )
        else:
            raise ValueError(
                f"Invalid optimization method: {RunConfig().optimization_method}"
                "Must be greedy or optimal."
            )
        # Try to fill coverage gaps with unassigned probes
        assigned = fill_coverage_gaps(self.df, assigned, ProbeConfig().spacing)

        if ProbeConfig().sequence_similarity > 0:
            match_percentage = ProbeConfig().sequence_similarity / 100
            assigned = self.filter_binding_probes(assigned, match_percentage)

        # Tm uniformity refinement: swap outlier probes for better alternatives
        assigned = self.refine_tm_uniformity(assigned, sequences)

        visualize_assignment(self.df, assigned, self.output()["coverage"].path)

        candidates = [sequence for sequence in sequences if sequence.id in assigned]
        util.log_and_check_candidates(
            self.logger, "OptimizeProbeCoverage", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["probes"].path, format="fasta")

    def _load_target_sequence(self) -> str | None:
        """Load the target gene sequence for accessibility scoring."""
        from .config import SequenceConfig

        sequence_fasta = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_sequence.fasta"
        )
        if not os.path.isfile(sequence_fasta):
            cfg = SequenceConfig()
            if cfg.sequence_file and os.path.isfile(cfg.sequence_file):
                sequence_fasta = cfg.sequence_file
            else:
                self.logger.debug("No target sequence for accessibility scoring")
                return None

        try:
            return str(next(Bio.SeqIO.parse(sequence_fasta, "fasta")).seq)
        except Exception:
            return None


def is_binding(seq1: str, seq2: str, match_percentage: float = 0.75) -> bool:
    """Check if seq1 is similar to (rev) complement of seq2 / if would bind.

    Uses k-mer overlap instead of full pairwise alignment for ~100x speedup.
    For probes of 20-25nt with k=8, this captures the biologically relevant
    contiguous complementarity that drives cross-hybridization.
    """
    if not (0 <= match_percentage <= 1):
        raise ValueError(
            "Matching percentage must be between 0 and 100. "
            f"Found `{match_percentage * 100}`"
        )
    min_len = min(len(seq1), len(seq2))
    if min_len < 3:
        return False

    # k scales with sequence length: k=8 for typical probes (20-25nt),
    # smaller for short sequences to maintain sensitivity
    k = min(8, max(3, min_len - 1))

    kmers1 = {seq1[i : i + k] for i in range(len(seq1) - k + 1)}
    n_kmers1 = len(kmers1)
    if n_kmers1 == 0:
        return False

    # Check complement and reverse complement
    comp = str(Bio.Seq.Seq(seq2).complement())
    rc = str(Bio.Seq.Seq(seq2).reverse_complement())

    for target in (comp, rc):
        kmers_target = {target[i : i + k] for i in range(len(target) - k + 1)}
        overlap = len(kmers1 & kmers_target) / max(n_kmers1, len(kmers_target))
        if overlap >= match_percentage:
            return True

    return False


def visualize_assignment(df: pd.DataFrame, assigned: list[str], filename: str) -> None:
    """Visualize the assignment as overlap matrix."""
    start = df["start"].min()
    end = df["end"].max()
    matrix = np.zeros((len(df) + 1, end - start))
    for idx, (_, row) in enumerate(df.iterrows()):
        matrix[idx + 1, row["start"] - start : row["end"] - start + 1] = (
            2 if row["name"] in assigned else 1
        )

    plt.figure()
    plt.title(f"Probe Coverage from {start} to {end}")
    plt.imshow(matrix)
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches="tight")
