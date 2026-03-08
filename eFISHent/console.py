"""Rich console utilities for eFISHent CLI."""

from contextlib import contextmanager
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from rich.theme import Theme

# Custom theme for eFISHent
THEME = Theme(
    {
        "info": "cyan",
        "warning": "yellow",
        "error": "red bold",
        "success": "green bold",
        "stage": "blue bold",
        "count": "magenta",
        "done": "green",
        "current": "cyan bold",
    }
)

# Global console instance
console = Console(theme=THEME)

# Module-level state for progress tracking
_pipeline_progress: Optional[Progress] = None
_pipeline_task_id: Optional[int] = None
_silent_mode: bool = False
_current_stage: int = 0
_total_stages: int = 8

# Filtering funnel data: list of (stage_name, candidate_count) tuples
_funnel_data: List[Tuple[str, int]] = []


def record_funnel_stage(stage_name: str, count: int) -> None:
    """Record a filtering stage result for the funnel visualization."""
    _funnel_data.append((stage_name, count))


def get_funnel_data() -> List[Tuple[str, int]]:
    """Get the accumulated funnel data."""
    return list(_funnel_data)


def reset_funnel_data() -> None:
    """Reset the funnel data (e.g. between runs)."""
    _funnel_data.clear()


def set_silent_mode(silent: bool) -> None:
    """Set silent mode to disable Rich output."""
    global _silent_mode
    _silent_mode = silent


def is_silent() -> bool:
    """Check if silent mode is enabled."""
    return _silent_mode


def get_rich_handler() -> RichHandler:
    """Get a configured Rich logging handler."""
    return RichHandler(
        console=console,
        show_time=True,
        show_path=False,
        rich_tracebacks=True,
        markup=True,
    )


def print_header(version: str) -> None:
    """Print the application header banner."""
    if _silent_mode:
        return
    console.print(
        Panel(
            f"[bold cyan]eFISHent v{version}[/bold cyan]",
            subtitle="RNA FISH probe designer",
            border_style="blue",
        )
    )


def print_completion(
    duration: str,
    output_files: Optional[List[str]] = None,
    summary: Optional[Dict] = None,
) -> None:
    """Print the completion message with design summary and output file locations.

    Args:
        duration: Human-readable duration string.
        output_files: List of output file paths.
        summary: Optional dict with keys: gene_name, probe_count, initial_count,
                 coverage_pct, length_range, length_median, tm_range, tm_median,
                 gc_range, gc_median.
    """
    if _silent_mode:
        return

    lines = []

    if summary:
        if summary.get("gene_name"):
            lines.append(f"  [bold]Gene:[/bold]      {summary['gene_name']}")
        if summary.get("probe_count") is not None:
            count_str = f"{summary['probe_count']}"
            if summary.get("initial_count"):
                count_str += f" selected from {summary['initial_count']:,} candidates"
            lines.append(f"  [bold]Probes:[/bold]    {count_str}")
        if summary.get("coverage_pct") is not None:
            lines.append(
                f"  [bold]Coverage:[/bold]  {summary['coverage_pct']:.1f}% of gene sequence"
            )
        if summary.get("length_range"):
            lo, hi = summary["length_range"]
            med = summary.get("length_median", "")
            med_str = f" (median {med})" if med else ""
            lines.append(f"  [bold]Length:[/bold]    {lo}-{hi} nt{med_str}")
        if summary.get("tm_range"):
            lo, hi = summary["tm_range"]
            med = summary.get("tm_median", "")
            med_str = f" (median {med:.1f}\u00b0C)" if med else ""
            lines.append(f"  [bold]TM range:[/bold]  {lo:.1f}-{hi:.1f}\u00b0C{med_str}")
        if summary.get("gc_range"):
            lo, hi = summary["gc_range"]
            med = summary.get("gc_median", "")
            med_str = f" (median {med:.1f}%)" if med else ""
            lines.append(f"  [bold]GC range:[/bold]  {lo:.1f}-{hi:.1f}%{med_str}")
        lines.append("")

    if output_files:
        lines.append("  [bold]Output:[/bold]")
        for path in output_files:
            lines.append(f"    [dim]\u2192[/dim] {path}")
        lines.append("")

    lines.append(f"  [dim]Completed in {duration}[/dim]")

    console.print(
        Panel(
            "\n".join(lines),
            title="[success]eFISHent \u2014 Design Complete[/success]",
            border_style="green",
        )
    )


def print_filtering_funnel() -> None:
    """Print a filtering funnel visualization showing probe counts at each stage."""
    if _silent_mode or not _funnel_data:
        return

    # Skip non-filtering stages (first = prepare, last = cleanup)
    # Keep: Generated, TM/GC, Alignment, K-mer, Secondary structure, Optimization
    stages = [
        (name, count)
        for name, count in _funnel_data
        if name not in ("Preparing gene sequence", "Finalizing output")
    ]
    if not stages:
        return

    max_count = max(count for _, count in stages)
    if max_count == 0:
        return

    # Short labels for compact display
    _short_names = {
        "Generating candidate probes": "Generated",
        "Filtering by TM/GC content": "TM/GC filter",
        "Aligning probes to genome": "Genome alignment",
        "Filtering by k-mer frequency": "K-mer filter",
        "Filtering by secondary structure": "Secondary structure",
        "Optimizing probe coverage": "Optimization",
    }

    bar_width = 24

    table = Table(
        show_header=False, box=None, padding=(0, 1), show_edge=False
    )
    table.add_column("Stage", style="dim", width=20, no_wrap=True)
    table.add_column("Bar", width=bar_width, no_wrap=True)
    table.add_column("Count", justify="right", width=7, no_wrap=True)
    table.add_column("Drop", width=10, no_wrap=True)

    for i, (stage_name, count) in enumerate(stages):
        filled = max(1, round(count / max_count * bar_width)) if max_count > 0 else 0
        bar = "\u2588" * filled
        short_name = _short_names.get(stage_name, stage_name)

        # Drop info
        if i == 0:
            drop_str = ""
        elif i == len(stages) - 1:
            drop_str = "[cyan]selected[/cyan]"
        else:
            prev_count = stages[i - 1][1]
            if prev_count > 0:
                drop_pct = (prev_count - count) / prev_count * 100
                drop_str = f"[dim]\u2193{drop_pct:>3.0f}%[/dim]"
            else:
                drop_str = ""

        # Color based on drop severity
        if i == 0:
            color = "green"
        elif i == len(stages) - 1:
            color = "cyan bold"
        else:
            prev_count = stages[i - 1][1]
            drop_pct = (prev_count - count) / prev_count * 100 if prev_count > 0 else 0
            if drop_pct > 50:
                color = "red"
            elif drop_pct > 20:
                color = "yellow"
            else:
                color = "green"

        table.add_row(
            short_name,
            f"[{color}]{bar}[/{color}]",
            f"{count:,}",
            drop_str,
        )

    console.print()
    console.print("[bold]Filtering Funnel[/bold]")
    console.print(table)


def print_stage(order: int, total: int, description: str) -> None:
    """Print a pipeline stage indicator with clear header."""
    global _current_stage, _total_stages
    if _silent_mode:
        return

    _current_stage = order
    _total_stages = total

    # Include candidate count from previous stage if available
    count_str = ""
    if _funnel_data:
        last_count = _funnel_data[-1][1]
        count_str = f" [dim]({last_count:,} probes)[/dim]"

    # Print a clear stage header
    console.print()
    console.print(
        f"[stage]Step {order}/{total}:[/stage] [current]{description}[/current]{count_str}"
    )

    # Update progress bar if active
    if _pipeline_progress is not None and _pipeline_task_id is not None:
        progress_desc = f"Step {order}/{total}: {description}"
        if _funnel_data:
            progress_desc += f" ({_funnel_data[-1][1]:,} probes)"
        _pipeline_progress.update(
            _pipeline_task_id,
            completed=order - 1,
            description=progress_desc,
        )


def print_candidate_count(name: str, count: int, count_prev: int = 0) -> None:
    """Print candidate count with optional comparison."""
    if _silent_mode:
        return
    if count_prev:
        filtered = count_prev - count
        console.print(
            f"  [done]\u2714[/done] [count]{count}[/count] candidates "
            f"[dim](filtered {filtered} from {count_prev})[/dim]"
        )
    else:
        console.print(f"  [done]\u2714[/done] [count]{count}[/count] candidates generated")

    # Update progress bar to show this stage is done
    if _pipeline_progress is not None and _pipeline_task_id is not None:
        _pipeline_progress.update(
            _pipeline_task_id,
            completed=_current_stage,
            description=f"Step {_current_stage}/{_total_stages} \u2714 {count:,} probes remaining",
        )


def print_parameter_warnings(warnings: List[str]) -> None:
    """Print parameter validation warnings in a panel."""
    if _silent_mode or not warnings:
        return

    lines = []
    for warning in warnings:
        # Indent continuation lines
        parts = warning.split("\n")
        lines.append(f"  [warning]\u26a0[/warning]  {parts[0]}")
        for part in parts[1:]:
            lines.append(f"    {part}")
        lines.append("")

    console.print(
        Panel(
            "\n".join(lines).rstrip(),
            title="[warning]Parameter Warnings[/warning]",
            border_style="yellow",
        )
    )


def print_warning(message: str) -> None:
    """Print a styled warning message."""
    if _silent_mode:
        return
    console.print(f"[warning]Warning:[/warning] {message}")


def print_error_panel(title: str, message: str, hint: Optional[str] = None) -> None:
    """Display an error in a formatted panel."""
    content = f"[error]{message}[/error]"
    if hint:
        content += f"\n\n[dim]Hint: {hint}[/dim]"
    console.print(Panel(content, title=f"[error]{title}[/error]", border_style="red"))


def _compute_probe_quality(row: pd.Series, df: pd.DataFrame) -> str:
    """Compute a quality score for a probe based on how central its metrics are.

    Returns a star rating string (3 stars max).
    """
    score = 3  # Start perfect, subtract for issues

    # TM centrality: penalize if near the edges of the range
    tm_range = df["TM"].max() - df["TM"].min()
    if tm_range > 0:
        tm_mid = (df["TM"].max() + df["TM"].min()) / 2
        tm_dist = abs(row["TM"] - tm_mid) / (tm_range / 2)
        if tm_dist > 0.8:
            score -= 1

    # GC balance: penalize extremes
    if row["GC"] < 30 or row["GC"] > 65:
        score -= 1

    # deltaG: penalize very negative (strong secondary structure)
    if row["deltaG"] < -8:
        score -= 1

    # Off-targets
    if "count" in row.index and row["count"] > 0:
        score -= 1

    score = max(1, min(3, score))
    return "\u2605" * score + "\u2606" * (3 - score)


def print_probe_table(df: pd.DataFrame, max_rows: int = 100) -> None:
    """Display probe data as a Rich table with quality scores."""
    if _silent_mode:
        return

    title = f"Probe Summary ({len(df)} probes)"

    table = Table(
        title=title, show_header=True, header_style="bold cyan", box=None
    )
    table.add_column("Name", style="green")
    table.add_column("Sequence", style="dim", max_width=25)
    table.add_column("Length", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("TM", justify="right")
    table.add_column("deltaG", justify="right")
    table.add_column("Score", justify="center")

    for _, row in df.head(max_rows).iterrows():
        seq = row["sequence"]
        seq_display = seq[:22] + "..." if len(seq) > 25 else seq
        quality = _compute_probe_quality(row, df)
        table.add_row(
            str(row["name"]),
            seq_display,
            str(row["length"]),
            f"{row['GC']:.1f}",
            f"{row['TM']:.1f}",
            f"{row['deltaG']:.1f}",
            quality,
        )

    if len(df) > max_rows:
        table.add_row("...", f"({len(df) - max_rows} more)", "", "", "", "", "")

    console.print(table)
    console.print()


@contextmanager
def pipeline_progress(total_stages: int = 8, mode: str = "pipeline"):
    """Context manager for pipeline progress tracking.

    Args:
        total_stages: Number of stages to track
        mode: "pipeline" for normal probe design, "analyze" for probe set analysis
    """
    global _pipeline_progress, _pipeline_task_id, _current_stage, _total_stages

    _current_stage = 0
    _total_stages = total_stages

    if _silent_mode:
        yield
        return

    console.print()
    if mode == "analyze":
        console.print("[dim]Analysis progress:[/dim]")
    else:
        console.print("[dim]Pipeline progress:[/dim]")

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=30),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
        transient=False,
    )

    with progress:
        _pipeline_progress = progress
        _pipeline_task_id = progress.add_task("Initializing...", total=total_stages)
        yield
        # Mark complete
        progress.update(
            _pipeline_task_id,
            completed=total_stages,
            description="[done]All steps completed![/done]",
        )

    _pipeline_progress = None
    _pipeline_task_id = None


def print_analysis_stage(step: int, total: int, description: str) -> None:
    """Print an analysis stage indicator."""
    global _current_stage, _total_stages
    if _silent_mode:
        return

    _current_stage = step
    _total_stages = total

    console.print(f"  [dim]\u2022[/dim] {description}")

    # Update progress bar if active
    if _pipeline_progress is not None and _pipeline_task_id is not None:
        _pipeline_progress.update(
            _pipeline_task_id,
            completed=step,
            description=f"Analyzing {description.lower()}...",
        )


def print_dependency_check(results: Dict[str, dict]) -> None:
    """Print dependency check results as a table."""
    table = Table(
        title="Dependency Check",
        show_header=True,
        header_style="bold cyan",
        border_style="dim",
    )
    table.add_column("Dependency", style="bold")
    table.add_column("Status", justify="center")
    table.add_column("Version / Path")
    table.add_column("Required For")

    all_ok = True
    for name, info in results.items():
        if info["found"]:
            status = "[green]\u2714[/green]"
            detail = info.get("version", info.get("path", ""))
        else:
            status = "[red]\u2718[/red]"
            detail = "[red]not found[/red]"
            all_ok = False
        table.add_row(name, status, detail, info.get("needed_for", ""))

    console.print(table)
    console.print()

    if all_ok:
        console.print("[success]All dependencies found![/success]")
    else:
        console.print(
            "[error]Missing dependencies detected.[/error] "
            "Run the install script for your platform:\n"
            "  [dim]macOS:[/dim]  ./install_macos.sh\n"
            "  [dim]Linux:[/dim]  ./install_linux.sh\n"
            "  [dim]Then:[/dim]   source ~/.local/efishent-deps/activate.sh"
        )
    console.print()


def print_missing_deps_error(missing: List[str]) -> None:
    """Print a concise error for missing dependencies during pipeline run."""
    deps = ", ".join(f"[bold]{d}[/bold]" for d in missing)
    console.print(
        Panel(
            f"[error]Missing required dependencies: {deps}[/error]\n\n"
            "Run the install script for your platform:\n"
            "  [dim]macOS:[/dim]  ./install_macos.sh\n"
            "  [dim]Linux:[/dim]  ./install_linux.sh\n"
            "  [dim]Then:[/dim]   source ~/.local/efishent-deps/activate.sh\n\n"
            "Or run [bold]efishent --check[/bold] for a full dependency report.",
            title="[error]Dependency Error[/error]",
            border_style="red",
        )
    )


@contextmanager
def spinner(description: str):
    """Context manager for spinner during subprocess operations."""
    if _silent_mode:
        yield
        return

    # If we're inside pipeline_progress, update the description to show subprocess
    if _pipeline_progress is not None and _pipeline_task_id is not None:
        old_desc = _pipeline_progress.tasks[_pipeline_task_id].description
        _pipeline_progress.update(_pipeline_task_id, description=f"[dim]{description}[/dim]")
        yield
        # Restore to show step progress
        _pipeline_progress.update(
            _pipeline_task_id,
            description=f"Running step {_current_stage}/{_total_stages}...",
        )
    else:
        # Standalone spinner - print inline status
        console.print(f"  [dim]\u23f3 {description}[/dim]", end="")
        yield
        console.print(" [done]\u2714[/done]")
