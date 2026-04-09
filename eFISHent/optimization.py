"""Find the optimal probeset with the highest coverage (fast or optimal)."""

from typing import List, Optional, Tuple
import functools
import logging
import multiprocessing
import os
import subprocess
import sys

import Bio.Seq
import Bio.SeqIO
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pe

from . import util
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .secondary_structure import SecondaryStructureFiltering


class OptimizeProbeCoverage(luigi.Task):
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

    def run_optimal_block(self, block: int) -> List[str]:
        """Run a single non-overlapping block."""
        df = self.df[self.df["block"] == block]
        return optimal_model(df, self.time_limit)

    def run_greedy(self):
        """Run sequential greedy model."""
        return greedy_model(self.df)

    def filter_binding_probes(
        self, assigned: List[str], match_percentage: float
    ) -> List[str]:
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
        self.df = util.create_data_table(sequences)
        self.df["end"] += ProbeConfig().spacing

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

    def _load_target_sequence(self) -> Optional[str]:
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

    def refine_tm_uniformity(
        self,
        assigned: List[str],
        sequences: List[Bio.SeqRecord.SeqRecord],
    ) -> List[str]:
        """Refine probe selection for Tm uniformity.

        For each probe whose Tm deviates > 5°C from the set median, check
        if an alternative probe at a nearby position has better Tm and swap.
        """
        from .basic_filtering import get_melting_temp

        if len(assigned) < 3:
            return assigned

        config = ProbeConfig()
        na = config.na_concentration
        formamide = config.formamide_concentration

        # Build sequence lookup
        seq_map = {s.id: s for s in sequences}

        # Compute Tm for all assigned probes
        assigned_tms = {}
        for name in assigned:
            if name in seq_map:
                assigned_tms[name] = get_melting_temp(
                    seq_map[name].seq, na, formamide
                )

        if not assigned_tms:
            return assigned

        median_tm = float(np.median(list(assigned_tms.values())))
        max_deviation = 5.0

        # Build position index of unassigned probes
        assigned_set = set(assigned)
        unassigned = [s for s in sequences if s.id not in assigned_set]

        swaps_made = 0
        new_assigned = list(assigned)

        for i, name in enumerate(new_assigned):
            if name not in assigned_tms:
                continue
            tm = assigned_tms[name]
            deviation = abs(tm - median_tm)
            if deviation <= max_deviation:
                continue

            # This probe is an outlier — look for better alternatives
            probe_row = self.df[self.df["name"] == name]
            if probe_row.empty:
                continue

            probe_start = int(name.split("-")[-1])
            probe_end = int(probe_row.iloc[0]["end"])

            best_alt, best_alt_deviation = _find_best_tm_alternative(
                unassigned, new_assigned, i, self.df,
                probe_start, probe_end, deviation,
                median_tm, na, formamide, get_melting_temp,
            )

            if best_alt and best_alt_deviation < deviation - 1.0:
                new_assigned[i] = best_alt.id
                swaps_made += 1

        if swaps_made > 0:
            self.logger.debug(
                f"Tm uniformity: swapped {swaps_made} probes for better Tm match "
                f"(median Tm: {median_tm:.1f}°C)"
            )

        return new_assigned


def _overlaps_other_assigned(
    alt_start: int,
    alt_end: int,
    new_assigned: List[str],
    skip_idx: int,
    df: pd.DataFrame,
) -> bool:
    """Check if an alternative probe overlaps any other assigned probe."""
    for j, other_name in enumerate(new_assigned):
        if j == skip_idx:
            continue
        other_row = df[df["name"] == other_name]
        if not other_row.empty:
            other_start = int(other_row.iloc[0]["start"])
            other_end = int(other_row.iloc[0]["end"])
            if is_overlapping((alt_start, alt_end), (other_start, other_end)):
                return True
    return False


def _find_best_tm_alternative(
    unassigned, new_assigned, skip_idx, df,
    probe_start, probe_end, current_deviation,
    median_tm, na, formamide, get_melting_temp_fn,
):
    """Find the best alternative probe with lower Tm deviation."""
    best_alt = None
    best_alt_deviation = current_deviation

    for alt in unassigned:
        alt_start = int(alt.id.split("-")[-1])
        alt_end = alt_start + len(alt)

        # Alternative must be at a nearby position (overlapping region)
        if alt_end < probe_start or alt_start > probe_end:
            continue

        if _overlaps_other_assigned(alt_start, alt_end, new_assigned, skip_idx, df):
            continue

        alt_tm = get_melting_temp_fn(alt.seq, na, formamide)
        alt_deviation = abs(alt_tm - median_tm)
        if alt_deviation < best_alt_deviation:
            best_alt = alt
            best_alt_deviation = alt_deviation

    return best_alt, best_alt_deviation


def is_overlapping(x: Tuple[int, int], y: Tuple[int, int]) -> bool:
    """Check if two ranges overlap. O(1) complexity."""
    return x[0] <= y[1] and y[0] <= x[1]


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


def fold_sequence(sequence: str) -> Optional[str]:
    """Fold a sequence and return the dot-bracket notation.

    Uses the bundled Fold binary (RNAstructure). Returns None if
    the binary is not available or folding fails.
    """
    from pathlib import Path

    file_path = Path(__file__).resolve().parent.as_posix()
    data_table = os.path.join(file_path, "data_tables/")
    if sys.platform in ("linux", "linux2"):
        fold_path = os.path.join(file_path, "Fold_linux")
    elif sys.platform == "darwin":
        fold_path = os.path.join(file_path, "Fold_osx")
    else:
        return None

    if not os.path.isfile(fold_path):
        return None

    os.environ["DATAPATH"] = data_table
    fasta_input = f">target\n{sequence}\n"
    try:
        result = subprocess.run(
            [fold_path, "-", "-", "--bracket", "--MFE"],
            input=fasta_input,
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 3:
            return lines[-1]
    except Exception:
        pass
    return None


def compute_accessibility_scores(
    df: pd.DataFrame, target_sequence: str, window: int = 200
) -> pd.Series:
    """Compute target RNA accessibility for each probe's binding site.

    Uses local RNA folding (200nt window) around each probe binding site
    to determine what fraction of the target is unpaired (accessible).
    Returns a Series of accessibility scores in [0, 1].
    """
    scores = []
    for _, row in df.iterrows():
        start = int(row["start"])
        end = int(row["end"])

        win_start = max(0, start - window // 2)
        win_end = min(len(target_sequence), end + window // 2)
        local_seq = target_sequence[win_start:win_end]

        if len(local_seq) < 10:
            scores.append(1.0)
            continue

        try:
            dot_bracket = fold_sequence(local_seq)
            if dot_bracket is None:
                scores.append(1.0)
                continue

            probe_start_in_window = start - win_start
            probe_end_in_window = end - win_start
            binding_region = dot_bracket[probe_start_in_window:probe_end_in_window]

            if len(binding_region) == 0:
                scores.append(1.0)
                continue

            unpaired = binding_region.count(".")
            accessibility = unpaired / len(binding_region)
            scores.append(round(accessibility, 4))
        except Exception:
            scores.append(1.0)

    return pd.Series(scores, index=df.index)


def compute_pre_optimization_quality(df: pd.DataFrame) -> pd.Series:
    """Compute a lightweight quality score for optimization weighting.

    Uses sequence properties available before the full quality pipeline:
    GC content, CpG fraction, and optionally accessibility.
    Data-driven from n=84 probe set analysis. Returns scores in [0, 1].
    """
    from .basic_filtering import get_gc_content, get_cpg_fraction

    # GC score: optimal 45-52%, steep penalty above 55%
    gc_vals = df["sequence"].apply(lambda s: get_gc_content(Bio.Seq.Seq(s)))

    def _gc_score(gc: float) -> float:
        if 45.0 <= gc <= 52.0:
            return 1.0
        elif gc < 45.0:
            return max(0.0, (gc - 20.0) / 25.0)
        elif gc <= 55.0:
            return 1.0 - (gc - 52.0) / 10.0
        else:
            return max(0.0, 0.7 - (gc - 55.0) / 15.0)

    gc_scores = gc_vals.apply(_gc_score)

    # CpG score: penalize high CpG fraction
    cpg_vals = df["sequence"].apply(lambda s: get_cpg_fraction(Bio.Seq.Seq(s)))
    cpg_scores = (1.0 - cpg_vals.clip(upper=0.10) / 0.10).clip(lower=0.0)

    # Accessibility score: prefer probes at unpaired (accessible) target sites
    has_accessibility = (
        "_accessibility" in df.columns
        and (df["_accessibility"] < 1.0).any()
    )

    if has_accessibility:
        acc_scores = df["_accessibility"].clip(lower=0.0, upper=1.0)
        # 50% GC, 30% CpG, 20% accessibility
        scores = (
            gc_scores * 0.50 + cpg_scores * 0.30 + acc_scores * 0.20
        ).clip(lower=0.1)
    else:
        # 60% GC, 40% CpG
        scores = (gc_scores * 0.6 + cpg_scores * 0.4).clip(lower=0.1)

    return scores


def greedy_model(df: pd.DataFrame) -> List[str]:
    """Run the greedy/fast model.

    Probes are sorted by start position, then by quality (highest first).
    At each position, the highest-quality non-overlapping probe is selected.
    """
    # Add quality scores for ranking
    df = df.copy()
    df["_quality"] = compute_pre_optimization_quality(df)

    df_sorted = df.sort_values(
        ["start", "_quality", "length"], ascending=[True, False, False]
    ).reset_index(drop=True)
    rec_array = df_sorted[["start", "end"]].to_records(index=False)
    probes = df_sorted["name"].values

    # Assign non-overlapping probes greedily
    assigned = []
    assigned_intervals = []

    for idx, probe in enumerate(probes):
        start, end = rec_array[idx]
        # Check if probe overlaps with ANY already assigned probe
        if all(
            not is_overlapping((start, end), interval)
            for interval in assigned_intervals
        ):
            assigned.append(probe)
            assigned_intervals.append((start, end))

    return assigned


def optimal_model(df: pd.DataFrame, time_limit: int) -> List[str]:
    """Run the quality-weighted optimal mathematical model.

    Maximizes coverage weighted by probe quality: max Σ quality[p] * covered[p,s].
    This prefers higher-quality probes when coverage is equal.
    """
    # Compute quality weights
    quality_weights = compute_pre_optimization_quality(df)
    quality_map = dict(zip(df["name"].values, quality_weights.values))

    # Convert dataframe to usable model inputs
    sequence = list(range(df["start"].min(), df["end"].max()))
    probes = df["name"].values
    probe_starts = {k: v for _, (k, v) in df[["name", "start"]].iterrows()}
    probe_ends = {k: v for _, (k, v) in df[["name", "end"]].iterrows()}

    coverages = {
        p: [t for t in sequence if t >= probe_starts[p] and t <= probe_ends[p]]
        for p in probes
    }

    # Model definition
    model = pe.ConcreteModel()
    model.sequence = pe.Set(initialize=sequence)
    model.probes = pe.Set(initialize=probes)
    model.covers = pe.Set(model.probes, within=model.sequence, initialize=coverages)
    model.covers_flat_set = pe.Set(
        initialize=[(p, s) for p in probes for s in model.covers[p]]
    )
    model.assign = pe.Var(model.probes, domain=pe.Binary)
    model.covered = pe.Var(model.covers_flat_set, domain=pe.Binary)

    # Quality-weighted objective: maximize Σ quality[p] * covered[p,s]
    obj = sum(
        quality_map.get(p, 1.0) * model.covered[p, s]
        for (p, s) in model.covers_flat_set
    )
    model.objective = pe.Objective(expr=obj, sense=pe.maximize)

    # Constraints
    def cover(model, p):
        return (
            sum(model.covered[p, s] for s in model.covers[p])
            == len(model.covers[p]) * model.assign[p]
        )

    def over_cover(model, s):
        cov_options = [(p, s) for p in model.probes if (p, s) in model.covers_flat_set]
        if not cov_options:
            return pe.Constraint.Skip
        return sum(model.covered[p, s] for (p, s) in cov_options) <= 1

    model.C1 = pe.Constraint(model.probes, rule=cover)
    model.C2 = pe.Constraint(model.sequence, rule=over_cover)
    solver = pe.SolverFactory("glpk")
    solver.options["tmlim"] = time_limit
    solver.solve(model, tee=False)

    assigned = [name for name, assign in model.assign.get_values().items() if assign]
    return assigned


def fill_coverage_gaps(
    df: pd.DataFrame, assigned: List[str], spacing: int
) -> List[str]:
    """Fill large coverage gaps with unassigned probes.

    After greedy/optimal assignment, there may be uncovered regions where no
    assigned probe lands. This function identifies gaps and tries to place
    the longest available unassigned probe in each gap.
    """
    if not assigned or len(df) == len(assigned):
        return assigned

    assigned_set = set(assigned)
    assigned_df = df[df["name"].isin(assigned_set)].sort_values("start")
    unassigned_df = df[~df["name"].isin(assigned_set)]

    if unassigned_df.empty:
        return assigned

    # Find gaps between consecutive assigned probes
    gaps = []
    prev_end = assigned_df.iloc[0]["end"]
    for _, row in assigned_df.iloc[1:].iterrows():
        gap_size = row["start"] - prev_end
        if gap_size > spacing:
            gaps.append((prev_end, row["start"]))
        prev_end = max(prev_end, row["end"])

    if not gaps:
        return assigned

    # For each gap, find the best unassigned probe that fits
    new_assigned = list(assigned)
    new_assigned_set = set(assigned)

    for gap_start, gap_end in gaps:
        # Find unassigned probes that fall within this gap
        candidates = unassigned_df[
            (unassigned_df["start"] >= gap_start)
            & (unassigned_df["end"] <= gap_end)
        ]
        if candidates.empty:
            continue

        # Pick the longest probe (most coverage)
        best = candidates.loc[candidates["length"].idxmax()]

        # Verify no overlap with any assigned probe
        overlaps = any(
            is_overlapping(
                (best["start"], best["end"]),
                (row["start"], row["end"]),
            )
            for _, row in df[df["name"].isin(new_assigned_set)].iterrows()
        )
        if not overlaps:
            new_assigned.append(best["name"])
            new_assigned_set.add(best["name"])

    return new_assigned


def visualize_assignment(df: pd.DataFrame, assigned: List[str], filename: str) -> None:
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
