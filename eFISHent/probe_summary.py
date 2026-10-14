"""Summarize the final probe set for the completion message."""

import os

import Bio.SeqIO
import pandas as pd

from . import util
from .config import SequenceConfig


class ProbeSummaryMixin:
    """Summary statistics (coverage, Tm/GC/length ranges) of the final probes."""

    def _get_gene_length(self) -> int:
        """Get the actual gene sequence length from the sequence file."""
        # Try the prepared sequence file (available before intermediate cleanup)
        sequence_fasta = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_sequence.fasta"
        )
        if os.path.isfile(sequence_fasta):
            seq = next(Bio.SeqIO.parse(sequence_fasta, "fasta"))
            return len(seq)
        # Try user-provided sequence file
        if SequenceConfig().sequence_file and os.path.isfile(SequenceConfig().sequence_file):
            seq = next(Bio.SeqIO.parse(SequenceConfig().sequence_file, "fasta"))
            return len(seq)
        return 0

    def _compute_summary(self, df: pd.DataFrame) -> dict:
        """Compute summary statistics from the final probe DataFrame."""
        from .console import get_funnel_data

        gene_name = util.get_gene_name(hashed=False)

        # Get actual gene length, fall back to probe span
        gene_length = self._get_gene_length()
        if gene_length == 0 and len(df) > 0:
            gene_length = int(df["end"].max())

        # Compute covered base pairs by merging overlapping intervals
        covered_bp = 0
        if len(df) > 0:
            intervals = sorted(zip(df["start"], df["end"]))
            merged_start, merged_end = intervals[0]
            for s, e in intervals[1:]:
                if s <= merged_end:
                    merged_end = max(merged_end, e)
                else:
                    covered_bp += merged_end - merged_start
                    merged_start, merged_end = s, e
            covered_bp += merged_end - merged_start

        coverage_pct = (covered_bp / gene_length * 100) if gene_length > 0 else 0.0

        # Get initial candidate count from funnel data
        funnel = get_funnel_data()
        initial_count = funnel[0][1] if funnel else None

        return {
            "gene_name": gene_name,
            "probe_count": len(df),
            "initial_count": initial_count,
            "coverage_pct": coverage_pct,
            "gene_length": gene_length,
            "length_range": (int(df["length"].min()), int(df["length"].max())),
            "length_median": int(df["length"].median()),
            "tm_range": (float(df["TM"].min()), float(df["TM"].max())),
            "tm_median": float(df["TM"].median()),
            "gc_range": (float(df["GC"].min()), float(df["GC"].max())),
            "gc_median": float(df["GC"].median()),
            "tm_std": float(df["TM"].std()) if len(df) > 1 else 0.0,
        }
