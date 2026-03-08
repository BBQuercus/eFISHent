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
    probe_df: Optional[pd.DataFrame] = None,
    verification: Optional[Dict] = None,
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
    if _silent_mode:
        return

    from rich.console import Group
    from rich.text import Text

    renderables = []

    # 1. Filtering funnel
    funnel_table = _build_funnel_table()
    if funnel_table is not None:
        renderables.append(Text(""))
        renderables.append(Text(" Filtering Funnel", style="bold"))
        renderables.append(funnel_table)

    # 2. Coverage map
    if probe_df is not None and summary and len(probe_df) > 0:
        gene_length = summary.get("gene_length", 0)
        coverage_pct = summary.get("coverage_pct", 0.0)
        coverage_map = _build_coverage_map(probe_df, gene_length, coverage_pct)
        if coverage_map is not None:
            renderables.append(Text(""))
            renderables.append(coverage_map)

    # 3. Probe table
    if probe_df is not None and len(probe_df) > 0:
        probe_table = _build_probe_table(probe_df, show_title=False)
        renderables.append(Text(""))
        renderables.append(probe_table)

    # 4. Summary stats
    if summary:
        renderables.append(Text(""))
        stats_table = Table(
            show_header=False, box=None, padding=(0, 1), show_edge=False
        )
        stats_table.add_column("Label", style="bold", no_wrap=True)
        stats_table.add_column("Value")

        if summary.get("gene_name"):
            stats_table.add_row("Gene:", summary["gene_name"])
        if summary.get("probe_count") is not None:
            count_str = f"{summary['probe_count']}"
            if summary.get("initial_count"):
                count_str += f" selected from {summary['initial_count']:,} candidates"
            stats_table.add_row("Probes:", count_str)
        if summary.get("coverage_pct") is not None:
            gene_len = summary.get("gene_length", 0)
            cov_str = f"{summary['coverage_pct']:.1f}%"
            if gene_len:
                cov_str += f" of {gene_len:,} nt"
            stats_table.add_row("Coverage:", cov_str)
        if summary.get("length_range"):
            lo, hi = summary["length_range"]
            med = summary.get("length_median", "")
            med_str = f" (median {med})" if med else ""
            stats_table.add_row("Length:", f"{lo}-{hi} nt{med_str}")
        if summary.get("tm_range"):
            lo, hi = summary["tm_range"]
            med = summary.get("tm_median", "")
            med_str = f" (median {med:.1f}\u00b0C)" if med else ""
            stats_table.add_row("TM range:", f"{lo:.1f}-{hi:.1f}\u00b0C{med_str}")
        if summary.get("gc_range"):
            lo, hi = summary["gc_range"]
            med = summary.get("gc_median", "")
            med_str = f" (median {med:.1f}%)" if med else ""
            stats_table.add_row("GC range:", f"{lo:.1f}-{hi:.1f}%{med_str}")
        renderables.append(stats_table)

    # 5. BLAST verification
    if verification:
        renderables.append(Text(""))
        renderables.append(_build_verification_summary(verification))

    # 6. Output files
    if output_files:
        renderables.append(Text(""))
        file_lines = [" [bold]Output:[/bold]"]
        for path in output_files:
            file_lines.append(f"   [dim]\u2192[/dim] {path}")
        renderables.append(Text.from_markup("\n".join(file_lines)))

    # 7. Duration
    renderables.append(Text(""))
    renderables.append(Text.from_markup(f" [dim]Completed in {duration}[/dim]"))

    console.print(
        Panel(
            Group(*renderables),
            title="[success]eFISHent \u2014 Design Complete[/success]",
            border_style="green",
        )
    )


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


def _build_funnel_table() -> Optional[Table]:
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


def _build_coverage_map(
    probe_df: pd.DataFrame,
    gene_length: int,
    coverage_pct: float,
) -> Optional["Text"]:
    """Build an ASCII coverage map showing probe binding positions on the gene.

    Returns a Rich Text renderable with a horizontal bar where filled blocks
    represent covered regions and light blocks represent gaps.
    """
    from rich.text import Text

    if gene_length <= 0:
        return None

    bar_width = 50

    # Build boolean coverage array across the gene
    coverage = [False] * gene_length
    for _, row in probe_df.iterrows():
        start = max(0, int(row["start"]))
        end = min(int(row["end"]), gene_length)
        for i in range(start, end):
            coverage[i] = True

    # Downsample to bar width — a cell is "covered" if any base in its range is
    bar = []
    for i in range(bar_width):
        region_start = int(i * gene_length / bar_width)
        region_end = int((i + 1) * gene_length / bar_width)
        bar.append(any(coverage[region_start:region_end]))

    # Build the visual bar
    text = Text()
    text.append(f" Coverage ({coverage_pct:.1f}% of {gene_length:,} nt)\n", style="bold")
    text.append("  5\u2032 ", style="dim")
    for covered in bar:
        if covered:
            text.append("\u2588", style="green")
        else:
            text.append("\u2591", style="dim")
    text.append(" 3\u2032\n", style="dim")

    # Scale line with position markers
    scale = [" "] * bar_width
    markers = [
        (0, "0"),
        (bar_width // 4, str(gene_length // 4)),
        (bar_width // 2, str(gene_length // 2)),
        (3 * bar_width // 4, str(3 * gene_length // 4)),
        (bar_width - 1, str(gene_length)),
    ]
    for pos, label in markers:
        # Place label at position, avoiding overlap
        for j, ch in enumerate(label):
            idx = pos + j
            if 0 <= idx < bar_width:
                scale[idx] = ch

    text.append("     " + "".join(scale).rstrip(), style="dim")
    return text


def _build_verification_summary(verification: Dict) -> "Text":
    """Build a Rich Text summary of BLAST verification results."""
    from rich.text import Text

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


def print_filtering_funnel() -> None:
    """Print a standalone filtering funnel visualization."""
    if _silent_mode or not _funnel_data:
        return
    table = _build_funnel_table()
    if table is None:
        return
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


def _build_probe_table(
    df: pd.DataFrame, max_rows: int = 100, show_title: bool = True
) -> Table:
    """Build a Rich Table for probe data. Returns the Table renderable."""
    title = f"Probe Summary ({len(df)} probes)" if show_title else None

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

    return table


def print_probe_table(df: pd.DataFrame, max_rows: int = 100) -> None:
    """Display probe data as a Rich table with quality scores."""
    if _silent_mode:
        return
    table = _build_probe_table(df, max_rows)
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
