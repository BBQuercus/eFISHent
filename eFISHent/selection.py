"""Select non-overlapping probes: weights, greedy/optimal models and gap filling."""

import Bio.Seq
import numpy as np
import pandas as pd
import pyomo.environ as pe


def is_overlapping(x: tuple[int, int], y: tuple[int, int]) -> bool:
    """Check if two ranges overlap. O(1) complexity."""
    return x[0] <= y[1] and y[0] <= x[1]


def compute_pre_optimization_quality(df: pd.DataFrame) -> pd.Series:
    """Compute a lightweight quality score for optimization weighting.

    Uses sequence properties available before the full quality pipeline:
    GC content, CpG fraction, and optionally accessibility.
    Data-driven from n=84 probe set analysis. Returns scores in [0, 1].
    """
    from .basic_filtering import get_cpg_fraction, get_gc_content

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


def cleanliness(df: pd.DataFrame) -> pd.Series:
    """Off-target cleanliness in [0, 1]: 1 for the lowest burden, 0 for the highest.

    Based on the burden rank among the candidates, because absolute burdens
    are dominated by chance 15-16 nt matches and differ a lot between genes.
    1 everywhere when no burden is known.
    """
    if "_burden" not in df or len(df) < 2 or df["_burden"].nunique() < 2:
        return pd.Series(1.0, index=df.index)
    rank = df["_burden"].rank(method="average")
    return 1.0 - (rank - rank.min()) / (rank.max() - rank.min())


def selection_weights(df: pd.DataFrame) -> pd.Series:
    """Per-probe weight for selection: probe count first, then specificity, then quality.

    Every probe is worth 1, plus up to 0.5 for off-target cleanliness and up
    to 0.25 for sequence quality (``_quality``). Because the bonus never
    exceeds 0.75, two probes always outweigh one, so the set is never shrunk
    for cleanliness; among arrangements with the same probe count the
    cleanest, best-quality one wins (TrueProbes-style: specific probes first).
    """
    quality = df["_quality"] if "_quality" in df else compute_pre_optimization_quality(df)
    return 1.0 + 0.5 * cleanliness(df) + 0.25 * quality


def greedy_model(df: pd.DataFrame) -> list[str]:
    """Select non-overlapping probes with the highest total weight.

    Exact weighted interval scheduling (dynamic programming over probes
    sorted by end, O(n log n)), using selection_weights(). Despite the name
    kept for the CLI option, the result is optimal for these weights; the
    previous left-to-right pass picked whichever probe started first.
    """
    if df.empty:
        return []
    df = df.copy()
    df["_quality"] = compute_pre_optimization_quality(df)
    weights = selection_weights(df).to_numpy(dtype=float)

    order = np.lexsort((df["start"].to_numpy(), df["end"].to_numpy()))
    starts = df["start"].to_numpy()[order]
    ends = df["end"].to_numpy()[order]
    weights = weights[order]
    names = df["name"].to_numpy()[order]

    # Probes overlap when start <= other end (is_overlapping), so probe i can
    # follow the first previous[i] probes in end order.
    previous = np.searchsorted(ends, starts, side="left")
    best = np.zeros(len(df) + 1)
    for i in range(len(df)):
        best[i + 1] = max(best[i], weights[i] + best[previous[i]])

    assigned = []
    i = len(df)
    while i > 0:
        if weights[i - 1] + best[previous[i - 1]] > best[i - 1]:
            assigned.append(names[i - 1])
            i = previous[i - 1]
        else:
            i -= 1
    return assigned[::-1]


def optimal_model(df: pd.DataFrame, time_limit: int) -> list[str]:
    """Run the quality-weighted optimal mathematical model.

    Maximizes coverage weighted by probe quality: max Σ quality[p] * covered[p,s].
    This prefers higher-quality probes when coverage is equal.
    """
    # Compute quality weights, discounted by off-target burden if known
    quality_weights = compute_pre_optimization_quality(df)
    quality_weights = quality_weights * (0.5 + 0.5 * cleanliness(df))
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
    df: pd.DataFrame, assigned: list[str], spacing: int
) -> list[str]:
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
