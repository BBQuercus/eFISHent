"""Probe set analysis: composite quality scores and the terminal summary."""

import pandas as pd

from .config import ProbeConfig
from .optimization import is_binding

KCAL_PER_MOL = "kcal/mol"

# Per-metric summary rows: (label, column, unit, decimals)
SUMMARY_METRICS = [
    ("Length", "length", "nt", 0),
    ("Tm", "TM", "°C", 1),
    ("GC", "GC", "%", 1),
    ("CpG fraction", "cpg_fraction", "", 3),
    ("ΔG (secondary)", "deltaG", KCAL_PER_MOL, 1),
    ("ΔG (binding)", "on_target_dg", KCAL_PER_MOL, 1),
    ("K-mer max", "kmers", "", 0),
    ("Homopolymer", "homopolymer_run", "nt", 0),
    ("Low complexity", "low_complexity", "", 3),
]


def compute_quality_scores(df: pd.DataFrame, cfg: ProbeConfig) -> pd.Series:
    """Compute composite quality scores matching the main pipeline."""
    scores = pd.DataFrame(index=df.index)

    # Tm score: distance from median (narrower spread = better)
    if len(df) > 1:
        median_tm = df["TM"].median()
        tm_range = (cfg.max_tm - cfg.min_tm) / 2 if cfg.max_tm > cfg.min_tm else 5.0
        scores["tm"] = 1.0 - (df["TM"] - median_tm).abs().clip(upper=tm_range) / tm_range
    else:
        scores["tm"] = 1.0

    # GC score
    def _gc_score(gc: float) -> float:
        if 45.0 <= gc <= 52.0:
            return 1.0
        elif gc < 45.0:
            return max(0.0, (gc - 20.0) / 25.0)
        elif gc <= 55.0:
            return 1.0 - (gc - 52.0) / 10.0
        else:
            return max(0.0, 0.7 - (gc - 55.0) / 15.0)

    scores["gc"] = df["GC"].apply(_gc_score)

    # CpG score
    scores["cpg"] = (1.0 - df["cpg_fraction"].clip(upper=0.10) / 0.10).clip(lower=0.0)

    # deltaG score (secondary structure)
    dg_max = abs(cfg.max_deltag) if cfg.max_deltag < 0 else 10.0
    scores["dg"] = 1.0 - df["deltaG"].abs().clip(upper=dg_max) / dg_max

    # K-mer score
    kmer_max = max(cfg.max_kmers, 1)
    scores["kmer"] = 1.0 - df["kmers"].clip(upper=kmer_max) / kmer_max

    # Off-target score
    if "off_targets" in df.columns:
        ot_max = max(cfg.max_off_targets + 1, 1)
        scores["ot"] = 1.0 - df["off_targets"].clip(upper=ot_max) / ot_max
    else:
        scores["ot"] = 1.0

    # Low-complexity score
    scores["lc"] = 1.0 - df["low_complexity"].clip(upper=0.3) / 0.3

    # Binding ΔG score
    dg_floor = -60.0
    scores["binding"] = (df["on_target_dg"].clip(upper=0.0) / dg_floor).clip(
        lower=0.0, upper=1.0
    )

    # Clip all to [0, 1]
    for col in scores.columns:
        scores[col] = scores[col].clip(lower=0.0, upper=1.0)

    # Weighted composite
    quality = (
        scores["tm"] * 0.18
        + scores["gc"] * 0.13
        + scores["cpg"] * 0.05
        + scores["dg"] * 0.13
        + scores["ot"] * 0.23
        + scores["kmer"] * 0.09
        + scores["lc"] * 0.09
        + scores["binding"] * 0.10
    )
    return (quality * 100).round(1)


def print_summary(df: pd.DataFrame, gene_len: int):
    """Print a text summary to the terminal."""
    from .console import console

    print_overview(df, gene_len)
    print_metrics(df)
    flags = collect_flags(df)

    if flags:
        console.print()
        console.print("  [bold yellow]Flags:[/bold yellow]")
        for f in flags:
            console.print(f"    [yellow]![/yellow] {f}")
    else:
        console.print()
        console.print("  [green]No issues detected[/green]")

    console.print()


def print_overview(df: pd.DataFrame, gene_len: int):
    """Print the probe count, gene length, and gene coverage."""
    from .console import console

    n = len(df)
    mapped = df[df["start"] >= 0]

    console.print()
    console.print(f"  [bold]Probe Set Analysis[/bold] — {n} probes")
    console.print(f"  Gene length: {gene_len:,} nt")
    if len(mapped) > 0:
        coverage = mapped.apply(
            lambda r: r["end"] - r["start"], axis=1
        ).sum() / gene_len * 100
        console.print(f"  Gene coverage: {coverage:.1f}%")
        console.print(
            f"  Probes mapped to gene: {len(mapped)}/{n}"
        )
    console.print()


def print_metrics(df: pd.DataFrame):
    """Print the per-metric median and range table."""
    from .console import console

    console.print("  [dim]Metric              Median     Range[/dim]")
    for label, col, unit, dec in SUMMARY_METRICS:
        if col not in df.columns:
            continue
        vals = df[col]
        fmt = f".{dec}f"
        med = f"{vals.median():{fmt}}"
        lo = f"{vals.min():{fmt}}"
        hi = f"{vals.max():{fmt}}"
        suffix = f" {unit}" if unit else ""
        console.print(
            f"  {label:<20s} {med:>8s}   {lo} – {hi}{suffix}"
        )

    if "off_targets" in df.columns:
        ot = df["off_targets"]
        console.print(
            f"  {'Off-targets':<20s} {ot.median():>8.0f}   {ot.min():.0f} – {ot.max():.0f}"
        )

    if "quality" in df.columns:
        q = df["quality"]
        console.print(
            f"  {'Quality score':<20s} {q.median():>8.1f}   {q.min():.1f} – {q.max():.1f}"
        )


def collect_flags(df: pd.DataFrame) -> list[str]:
    """Collect warnings about potentially problematic probes."""
    flags = []
    g4_count = df["g_quadruplex"].sum()
    if g4_count > 0:
        flags.append(f"{g4_count} probe(s) contain G-quadruplex motifs")
    high_gc = (df["GC"] > 55).sum()
    if high_gc > 0:
        flags.append(f"{high_gc} probe(s) have GC > 55% (risk of non-specific binding)")
    high_lc = (df["low_complexity"] > 0.1).sum()
    if high_lc > 0:
        flags.append(f"{high_lc} probe(s) have elevated low-complexity score")
    if "off_targets" in df.columns:
        multi_ot = (df["off_targets"] > 0).sum()
        if multi_ot > 0:
            flags.append(f"{multi_ot} probe(s) have genome off-target hits")

    # Cross-hybridization
    binding_pairs = count_binding_pairs(df["sequence"].values)
    if binding_pairs > 0:
        flags.append(f"{binding_pairs} probe pair(s) may cross-hybridize (≥75% k-mer overlap)")
    return flags


def count_binding_pairs(seqs) -> int:
    """Count probe pairs that may cross-hybridize (≥75% k-mer overlap)."""
    binding_pairs = 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            if is_binding(seqs[i], seqs[j], 0.75):
                binding_pairs += 1
    return binding_pairs
