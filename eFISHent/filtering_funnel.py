"""Filtering funnel: per-stage candidate counts and their Rich visualization."""

from rich.table import Table

from .console_base import console, is_silent

# Filtering funnel data: list of (stage_name, candidate_count) tuples
_funnel_data: list[tuple[str, int]] = []


def record_funnel_stage(stage_name: str, count: int) -> None:
    """Record a filtering stage result for the funnel visualization."""
    _funnel_data.append((stage_name, count))


def get_funnel_data() -> list[tuple[str, int]]:
    """Get the accumulated funnel data."""
    return list(_funnel_data)


def reset_funnel_data() -> None:
    """Reset the funnel data (e.g. between runs)."""
    _funnel_data.clear()


# Short labels for funnel display
_FUNNEL_SHORT_NAMES = {
    "Generating candidate probes": "Generated",
    "Filtering by TM/GC content": "TM/GC filter",
    "Aligning probes to genome": "Genome alignment",
    "Filtering by k-mer frequency": "K-mer filter",
    "Filtering by secondary structure": "Secondary structure",
    "Optimizing probe coverage": "Optimization",
}

# Stages to skip in the funnel (non-filtering)
_FUNNEL_SKIP_STAGES = {"Preparing gene sequence", "Finalizing output"}


def _get_drop_pct(count: int, prev_count: int) -> float:
    """Calculate percentage drop between two stage counts."""
    return (prev_count - count) / prev_count * 100 if prev_count > 0 else 0


def _get_stage_color(drop_pct: float) -> str:
    """Return Rich color based on drop severity."""
    if drop_pct > 50:
        return "red"
    if drop_pct > 20:
        return "yellow"
    return "green"


def _build_funnel_table() -> Table | None:
    """Build a Rich Table for the filtering funnel. Returns None if no data."""
    stages = [
        (name, count)
        for name, count in _funnel_data
        if name not in _FUNNEL_SKIP_STAGES
    ]
    if not stages:
        return None

    max_count = max(count for _, count in stages)
    if max_count == 0:
        return None

    bar_width = 24
    last = len(stages) - 1

    table = Table(
        show_header=False, box=None, padding=(0, 1), show_edge=False,
        expand=False,
    )
    table.add_column("Stage", style="dim", min_width=20, no_wrap=True)
    table.add_column("Bar", min_width=bar_width, no_wrap=True)
    table.add_column("Count", justify="right", min_width=7, no_wrap=True)
    table.add_column("Drop", min_width=8, no_wrap=True)

    for i, (stage_name, count) in enumerate(stages):
        filled = max(1, round(count / max_count * bar_width))
        bar = "\u2588" * filled
        short_name = _FUNNEL_SHORT_NAMES.get(stage_name, stage_name)
        drop_pct = _get_drop_pct(count, stages[i - 1][1]) if i > 0 else 0

        if i == 0:
            color, drop_str = "green", ""
        elif i == last:
            color, drop_str = "cyan bold", "[cyan]selected[/cyan]"
        else:
            color = _get_stage_color(drop_pct)
            drop_str = f"[dim]\u2193{drop_pct:>3.0f}%[/dim]" if drop_pct else ""

        table.add_row(
            short_name,
            f"[{color}]{bar}[/{color}]",
            f"{count:,}",
            drop_str,
        )

    return table


def print_filtering_funnel() -> None:
    """Print a standalone filtering funnel visualization."""
    if is_silent() or not _funnel_data:
        return
    table = _build_funnel_table()
    if table is None:
        return
    console.print()
    console.print("[bold]Filtering Funnel[/bold]")
    console.print(table)
