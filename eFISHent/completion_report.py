"""Completion report printed after a successful probe design run."""

import pandas as pd
from rich.table import Table
from rich.text import Text

from .console_base import console, is_silent
from .filtering_funnel import _build_funnel_table
from .probe_display import _build_coverage_map, _build_probe_table


def print_completion(
    duration: str,
    output_files: list[str] | None = None,
    summary: dict | None = None,
    probe_df: pd.DataFrame | None = None,
    verification: dict | None = None,
) -> None:
    """Print the completion panel with funnel, coverage map, probe table, and summary.

    Args:
        duration: Human-readable duration string.
        output_files: List of output file paths.
        summary: Optional dict with keys: gene_name, probe_count, initial_count,
                 coverage_pct, gene_length, length_range, length_median,
                 tm_range, tm_median, gc_range, gc_median.
        probe_df: Optional DataFrame with probe data to show in the table.
        verification: Optional dict with BLAST verification results.
    """
    if is_silent():
        return

    from rich.rule import Rule

    console.print()
    console.print(Rule("[success]Design Complete[/success]", style="green"))

    # 1. Filtering funnel
    _print_funnel_section()

    has_probes = probe_df is not None and len(probe_df) > 0

    # 2. Coverage map
    if has_probes and summary:
        _print_coverage_section(probe_df, summary)

    # 3. Probe table
    if has_probes:
        console.print()
        console.print(_build_probe_table(probe_df, show_title=False))

    # 4. Summary stats
    if summary:
        console.print()
        console.print(_build_stats_table(summary, probe_df))

    # 5. BLAST verification
    if verification:
        console.print()
        console.print(_build_verification_summary(verification))

    # 6. Output files
    if output_files:
        _print_output_files(output_files)

    # 7. Duration + closing rule
    console.print()
    console.print(f"  [dim]Completed in {duration}[/dim]")
    console.print(Rule(style="green"))


def _print_funnel_section() -> None:
    """Print the filtering funnel section of the completion panel, if any."""
    funnel_table = _build_funnel_table()
    if funnel_table is None:
        return
    console.print()
    console.print("  [bold]Filtering Funnel[/bold]")
    console.print(funnel_table)


def _print_coverage_section(probe_df: pd.DataFrame, summary: dict) -> None:
    """Print the coverage map section of the completion panel, if any."""
    gene_length = summary.get("gene_length", 0)
    coverage_pct = summary.get("coverage_pct", 0.0)
    coverage_map = _build_coverage_map(probe_df, gene_length, coverage_pct)
    if coverage_map is None:
        return
    console.print()
    console.print(coverage_map)


def _print_output_files(output_files: list[str]) -> None:
    """Print the list of output files."""
    console.print()
    console.print("  [bold]Output[/bold]")
    for path in output_files:
        console.print(f"    [dim]\u2192[/dim] {path}")


def _format_gene_stat(summary: dict) -> str | None:
    """Format the gene name summary row."""
    return summary.get("gene_name") or None


def _format_probe_count_stat(summary: dict) -> str | None:
    """Format the probe count summary row."""
    if summary.get("probe_count") is None:
        return None
    count_str = f"{summary['probe_count']}"
    if summary.get("initial_count"):
        count_str += f" selected from {summary['initial_count']:,} candidates"
    return count_str


def _format_coverage_stat(summary: dict) -> str | None:
    """Format the coverage summary row."""
    if summary.get("coverage_pct") is None:
        return None
    gene_len = summary.get("gene_length", 0)
    cov_str = f"{summary['coverage_pct']:.1f}%"
    if gene_len:
        cov_str += f" of {gene_len:,} nt"
    return cov_str


def _format_length_stat(summary: dict) -> str | None:
    """Format the probe length summary row."""
    if not summary.get("length_range"):
        return None
    lo, hi = summary["length_range"]
    med = summary.get("length_median", "")
    med_str = f" (median {med})" if med else ""
    return f"{lo}-{hi} nt{med_str}"


def _format_tm_stat(summary: dict) -> str | None:
    """Format the melting temperature summary row."""
    if not summary.get("tm_range"):
        return None
    lo, hi = summary["tm_range"]
    med = summary.get("tm_median", "")
    med_str = f" (median {med:.1f}\u00b0C)" if med else ""
    tm_std = summary.get("tm_std")
    if tm_std is not None:
        med_str += f", \u03c3 {tm_std:.1f}\u00b0C"
    return f"{lo:.1f}-{hi:.1f}\u00b0C{med_str}"


def _format_gc_stat(summary: dict) -> str | None:
    """Format the GC content summary row."""
    if not summary.get("gc_range"):
        return None
    lo, hi = summary["gc_range"]
    med = summary.get("gc_median", "")
    med_str = f" (median {med:.1f}%)" if med else ""
    return f"{lo:.1f}-{hi:.1f}%{med_str}"


# Summary stat rows: (label, formatter returning the value or None to skip)
_SUMMARY_STAT_ROWS = (
    ("Gene:", _format_gene_stat),
    ("Probes:", _format_probe_count_stat),
    ("Coverage:", _format_coverage_stat),
    ("Length:", _format_length_stat),
    ("TM range:", _format_tm_stat),
    ("GC range:", _format_gc_stat),
)

_RECOMMENDATION_STYLES = {"PASS": "green", "FLAG": "yellow", "FAIL": "red"}


def _format_recommendation_breakdown(probe_df: pd.DataFrame) -> str | None:
    """Summarize PASS/FLAG/FAIL recommendation counts from probe data."""
    if "recommendation" not in probe_df.columns:
        return None
    counts = probe_df["recommendation"].value_counts()
    parts = []
    for label, style in _RECOMMENDATION_STYLES.items():
        n = counts.get(label, 0)
        if n > 0:
            parts.append(f"[{style}]{n} {label}[/{style}]")
    return ", ".join(parts) if parts else None


def _format_off_target_summary(probe_df: pd.DataFrame) -> str | None:
    """Summarize transcriptome off-targets from probe data."""
    if "txome_off_targets" not in probe_df.columns:
        return None
    total_ot = int(probe_df["txome_off_targets"].sum())
    probes_with_ot = int((probe_df["txome_off_targets"] > 0).sum())
    if probes_with_ot > 0:
        return f"{probes_with_ot}/{len(probe_df)} probes have transcriptome off-targets ({total_ot} total)"
    return "[green]none detected[/green]"


def _build_stats_table(summary: dict, probe_df: pd.DataFrame | None) -> Table:
    """Build the summary statistics table of the completion panel."""
    stats_table = Table(
        show_header=False, box=None, padding=(0, 1), show_edge=False
    )
    stats_table.add_column("Label", style="bold", no_wrap=True)
    stats_table.add_column("Value")

    rows = [(label, formatter(summary)) for label, formatter in _SUMMARY_STAT_ROWS]
    if probe_df is not None:
        # Recommendation breakdown and off-target summary from probe data
        rows.append(("Quality:", _format_recommendation_breakdown(probe_df)))
        rows.append(("Off-targets:", _format_off_target_summary(probe_df)))

    for label, value in rows:
        if value is not None:
            stats_table.add_row(label, value)
    return stats_table


def _build_verification_summary(verification: dict) -> Text:
    """Build a Rich Text summary of BLAST verification results."""
    total = verification["total"]
    clean = verification["clean"]
    flagged = verification.get("flagged", {})
    max_show = 5

    text = Text()
    text.append(" BLAST Verification\n", style="bold")

    if clean == total:
        text.append("  \u2714 ", style="green")
        text.append(f"{clean}/{total} probes verified ", style="green")
        text.append("(independent BLAST cross-check)", style="dim")
    else:
        text.append("  \u2714 ", style="green")
        text.append(f"{clean}/{total} probes clean", style="green")
        n_flagged = total - clean
        text.append(f", {n_flagged} with additional BLAST hits:\n", style="yellow")
        max_expected = verification.get("max_expected", 1)
        items = sorted(flagged.items(), key=lambda x: -x[1])
        for probe_name, hit_count in items[:max_show]:
            text.append(f"    \u26a0 {probe_name}: ", style="yellow")
            text.append(f"{hit_count} loci ", style="yellow")
            text.append(f"(expected \u2264{max_expected})\n", style="dim")
        if len(items) > max_show:
            text.append(
                f"    [dim]...and {len(items) - max_show} more[/dim]",
            )

    return text
