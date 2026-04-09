"""Rich console utilities for eFISHent CLI."""

from contextlib import contextmanager
from typing import Dict, List, Optional, Tuple
import logging
import time

import pandas as pd
from rich.console import Console, Group
from rich.live import Live
from rich.text import Text
from rich.panel import Panel
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
_silent_mode: bool = False
_current_stage: int = 0
_total_stages: int = 8

# Live display state
_live: Optional[Live] = None
_steps: List[Dict] = []  # Each: {name, status, result, elapsed, start_time}
_step_status: str = ""  # Transient status line for current step
_pipeline_mode: str = "pipeline"

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


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds into a compact string."""
    if seconds < 1:
        return "<1s"
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    return f"{minutes}m {secs:02d}s"


# Braille spinner frames
_SPINNER_FRAMES = "\u280b\u2819\u2839\u2838\u283c\u2834\u2826\u2827\u2807\u280f"


class _PipelineDisplay:
    """Dynamic renderable that rebuilds on every Rich render cycle.

    This ensures elapsed times tick and the spinner animates without
    needing explicit refresh calls.
    """

    def __init__(self) -> None:
        self._frame: int = 0

    def __rich_console__(self, console: Console, options: object):  # type: ignore[override]
        self._frame += 1
        spinner_char = _SPINNER_FRAMES[self._frame % len(_SPINNER_FRAMES)]

        table = Table(
            show_header=False,
            box=None,
            padding=(0, 1),
            show_edge=False,
            expand=True,
        )
        table.add_column("Icon", width=2, no_wrap=True)
        table.add_column("Step", no_wrap=True)
        table.add_column("Spacer", ratio=1)  # Absorbs remaining width
        table.add_column("Result", no_wrap=True, justify="right")
        table.add_column("Time", justify="right", width=8, no_wrap=True)

        has_running = False
        for step in _steps:
            status = step["status"]
            name = step["name"]
            result = step.get("result", "")

            if status == "done":
                icon = "[green]\u2714[/green]"
                elapsed = step.get("elapsed", "")
                name_style = "dim"
                result_style = "dim"
            elif status == "running":
                has_running = True
                icon = f"[cyan]{spinner_char}[/cyan]"
                # Compute elapsed live
                elapsed = _format_elapsed(time.time() - step["start_time"])
                name_style = "bold"
                result_style = "cyan"
            else:
                icon = " "
                elapsed = ""
                name_style = "dim"
                result_style = "dim"

            table.add_row(
                icon,
                f"[{name_style}]{name}[/{name_style}]",
                "",
                f"[{result_style}]{result}[/{result_style}]" if result else "",
                f"[dim]{elapsed}[/dim]" if elapsed else "",
            )

            # Show transient status under the running step
            if status == "running" and _step_status:
                table.add_row(
                    "",
                    f"  [cyan]{spinner_char}[/cyan] [dim]{_step_status}[/dim]",
                    "",
                    "",
                    "",
                )

        # Show status line when no step is running (e.g. dependency tasks)
        if not has_running and _step_status:
            table.add_row(
                f"[cyan]{spinner_char}[/cyan]",
                f"[dim]{_step_status}[/dim]",
                "",
                "",
                "",
            )

        # Total elapsed at the bottom
        if _steps:
            first_start = _steps[0].get("start_time")
            if first_start is not None:
                total = _format_elapsed(time.time() - first_start)
                table.add_row("", "", "", "", "")
                table.add_row(
                    "",
                    "",
                    "",
                    f"[dim]{total} elapsed[/dim]",
                    "",
                )

        yield table


# Singleton display instance — reused across pipeline runs
_pipeline_display = _PipelineDisplay()


def _refresh_live() -> None:
    """Force a refresh of the live display if active."""
    if _live is not None:
        _live.refresh()


class LiveLogHandler(logging.Handler):
    """Logging handler that routes messages to the live step display."""

    def emit(self, record: logging.LogRecord) -> None:
        global _step_status
        if _live is not None and not _silent_mode:
            # Show the message as transient status under current step
            msg = record.getMessage()
            # Strip leading whitespace for indented messages (like transcript list)
            msg = msg.strip()
            if msg:
                _step_status = msg
                _refresh_live()
        elif not _silent_mode:
            # Outside live display, print directly
            console.print(f"  [dim]{record.getMessage()}[/dim]")


def get_live_log_handler() -> LiveLogHandler:
    """Get a logging handler that routes to the live display."""
    return LiveLogHandler()


def get_rich_handler() -> logging.Handler:
    """Get a configured logging handler.

    Returns a LiveLogHandler that integrates with the live step display.
    """
    return LiveLogHandler()


def print_header(version: str) -> None:
    """Print the application header banner."""
    if _silent_mode:
        return
    from rich.rule import Rule

    console.print()
    console.print(f"  [bold cyan]eFISHent v{version}[/bold cyan]  [dim]RNA FISH probe designer[/dim]")
    console.print(Rule(style="blue"))


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

    renderables = []

    from rich.rule import Rule

    console.print()
    console.print(Rule("[success]Design Complete[/success]", style="green"))

    # 1. Filtering funnel
    funnel_table = _build_funnel_table()
    if funnel_table is not None:
        console.print()
        console.print("  [bold]Filtering Funnel[/bold]")
        console.print(funnel_table)

    # 2. Coverage map
    if probe_df is not None and summary and len(probe_df) > 0:
        gene_length = summary.get("gene_length", 0)
        coverage_pct = summary.get("coverage_pct", 0.0)
        coverage_map = _build_coverage_map(probe_df, gene_length, coverage_pct)
        if coverage_map is not None:
            console.print()
            console.print(coverage_map)

    # 3. Probe table
    if probe_df is not None and len(probe_df) > 0:
        probe_table = _build_probe_table(probe_df, show_title=False)
        console.print()
        console.print(probe_table)

    # 4. Summary stats
    if summary:
        console.print()
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
            tm_std = summary.get("tm_std")
            if tm_std is not None:
                med_str += f", \u03c3 {tm_std:.1f}\u00b0C"
            stats_table.add_row("TM range:", f"{lo:.1f}-{hi:.1f}\u00b0C{med_str}")
        if summary.get("gc_range"):
            lo, hi = summary["gc_range"]
            med = summary.get("gc_median", "")
            med_str = f" (median {med:.1f}%)" if med else ""
            stats_table.add_row("GC range:", f"{lo:.1f}-{hi:.1f}%{med_str}")

        # Recommendation breakdown from probe data
        if probe_df is not None and "recommendation" in probe_df.columns:
            counts = probe_df["recommendation"].value_counts()
            parts = []
            for label in ("PASS", "FLAG", "FAIL"):
                n = counts.get(label, 0)
                if n > 0:
                    style = {"PASS": "green", "FLAG": "yellow", "FAIL": "red"}[label]
                    parts.append(f"[{style}]{n} {label}[/{style}]")
            if parts:
                stats_table.add_row("Quality:", ", ".join(parts))

        # Off-target summary from probe data
        if probe_df is not None and "txome_off_targets" in probe_df.columns:
            total_ot = int(probe_df["txome_off_targets"].sum())
            probes_with_ot = int((probe_df["txome_off_targets"] > 0).sum())
            if probes_with_ot > 0:
                stats_table.add_row(
                    "Off-targets:",
                    f"{probes_with_ot}/{len(probe_df)} probes have transcriptome off-targets ({total_ot} total)",
                )
            else:
                stats_table.add_row("Off-targets:", "[green]none detected[/green]")

        console.print(stats_table)

    # 5. BLAST verification
    if verification:
        console.print()
        console.print(_build_verification_summary(verification))

    # 6. Output files
    if output_files:
        console.print()
        console.print("  [bold]Output[/bold]")
        for path in output_files:
            console.print(f"    [dim]\u2192[/dim] {path}")

    # 7. Duration + closing rule
    console.print()
    console.print(f"  [dim]Completed in {duration}[/dim]")
    console.print(Rule(style="green"))


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
) -> Optional[Text]:
    """Build an ASCII coverage map showing probe binding positions on the gene.

    Returns a Rich Text renderable with a horizontal bar where filled blocks
    represent covered regions and light blocks represent gaps.
    """
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

    # Scale line with position markers at 0%, 25%, 50%, 75%, 100%
    # Use a wider buffer to avoid truncation at the right edge
    scale_width = bar_width + 10
    scale = [" "] * scale_width
    markers = [
        (0, "0"),
        (bar_width // 4, str(gene_length // 4)),
        (bar_width // 2, str(gene_length // 2)),
        (3 * bar_width // 4, str(3 * gene_length // 4)),
    ]
    # Right-align the last label so it doesn't overflow
    last_label = str(gene_length)
    markers.append((bar_width - len(last_label), last_label))

    for pos, label in markers:
        for j, ch in enumerate(label):
            idx = pos + j
            if 0 <= idx < scale_width:
                scale[idx] = ch

    text.append("     " + "".join(scale).rstrip(), style="dim")
    return text


def _build_verification_summary(verification: Dict) -> Text:
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


def add_cached_step(order: int, total: int, description: str) -> None:
    """Add a step that was already completed (cached by Luigi)."""
    if _silent_mode:
        return

    if order > 0:
        step_name = f"Step {order}/{total}: {description}"
    else:
        step_name = description

    _steps.append({
        "name": step_name,
        "status": "done",
        "result": "[dim]cached[/dim]",
        "elapsed": "",
        "start_time": 0,
    })
    _refresh_live()


def print_stage(order: int, total: int, description: str) -> None:
    """Print a pipeline stage indicator — updates the live step display."""
    global _current_stage, _total_stages, _step_status
    if _silent_mode:
        return

    _current_stage = order
    _total_stages = total

    # Finish the previous running step
    for step in _steps:
        if step["status"] == "running":
            step["status"] = "done"
            step["elapsed"] = _format_elapsed(time.time() - step["start_time"])

    # Include candidate count from previous stage if available
    count_str = ""
    if _funnel_data:
        last_count = _funnel_data[-1][1]
        count_str = f"({last_count:,} probes)"

    # Add new step
    if order > 0:
        step_name = f"Step {order}/{total}: {description}"
    else:
        step_name = description

    _steps.append({
        "name": step_name,
        "status": "running",
        "result": count_str,
        "elapsed": "",
        "start_time": time.time(),
    })
    _step_status = ""
    _refresh_live()


def print_candidate_count(_name: str, count: int, count_prev: int = 0) -> None:
    """Print candidate count — updates the current step's result in the live display."""
    global _step_status
    if _silent_mode:
        return

    # Update the current running step's result
    for step in reversed(_steps):
        if step["status"] == "running":
            if count_prev:
                filtered = count_prev - count
                drop_pct = (filtered / count_prev * 100) if count_prev else 0
                step["result"] = f"{count:,} remaining  (\u2212{drop_pct:.1f}%)"
            else:
                step["result"] = f"{count:,} probes"
            break

    _step_status = ""
    _refresh_live()


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
    # If live display is active, show as transient status
    if _live is not None:
        global _step_status
        _step_status = f"\u26a0 {message}"
        _refresh_live()
    else:
        console.print(f"[warning]Warning:[/warning] {message}")


def print_error_panel(title: str, message: str, hint: Optional[str] = None) -> None:
    """Display an error in a formatted panel."""
    # Stop live display before showing error so it doesn't interfere
    if _live is not None:
        _live.stop()
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
    has_recommendation = "recommendation" in df.columns
    has_off_targets = "txome_off_targets" in df.columns

    table.add_column("Name", style="green")
    table.add_column("Sequence", style="dim", max_width=25)
    table.add_column("Length", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("TM", justify="right")
    table.add_column("deltaG", justify="right")
    if has_off_targets:
        table.add_column("Off-targets", justify="right")
    if has_recommendation:
        table.add_column("Rec.", justify="center")
    table.add_column("Score", justify="center")

    for _, row in df.head(max_rows).iterrows():
        seq = row["sequence"]
        seq_display = seq[:22] + "..." if len(seq) > 25 else seq
        quality = _compute_probe_quality(row, df)
        cells = [
            str(row["name"]),
            seq_display,
            str(row["length"]),
            f"{row['GC']:.1f}",
            f"{row['TM']:.1f}",
            f"{row['deltaG']:.1f}",
        ]
        if has_off_targets:
            ot = int(row.get("txome_off_targets", 0))
            cells.append(str(ot) if ot > 0 else "[dim]0[/dim]")
        if has_recommendation:
            rec = str(row.get("recommendation", ""))
            rec_style = {"PASS": "[green]PASS[/green]", "FLAG": "[yellow]FLAG[/yellow]", "FAIL": "[red]FAIL[/red]"}
            cells.append(rec_style.get(rec, rec))
        cells.append(quality)
        table.add_row(*cells)

    if len(df) > max_rows:
        ncols = 7 + int(has_off_targets) + int(has_recommendation)
        filler = [""] * (ncols - 2)
        table.add_row("...", f"({len(df) - max_rows} more)", *filler)

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
    """Context manager for pipeline progress tracking using a live step display.

    Args:
        total_stages: Number of stages to track
        mode: "pipeline" for normal probe design, "analyze" for probe set analysis
    """
    global _live, _steps, _step_status, _current_stage, _total_stages, _pipeline_mode

    _current_stage = 0
    _total_stages = total_stages
    _steps = []
    _step_status = ""
    _pipeline_mode = mode

    if _silent_mode:
        yield
        return

    live = Live(
        _pipeline_display,
        console=console,
        refresh_per_second=8,
        transient=False,
    )

    with live:
        _live = live
        yield
        # Finish the last running step
        for step in _steps:
            if step["status"] == "running":
                step["status"] = "done"
                step["elapsed"] = _format_elapsed(time.time() - step["start_time"])
        _refresh_live()

    _live = None


def print_analysis_stage(step: int, total: int, description: str) -> None:
    """Print an analysis stage indicator — uses the live step display."""
    global _current_stage, _total_stages, _step_status
    if _silent_mode:
        return

    _current_stage = step
    _total_stages = total

    # Finish the previous running step
    for s in _steps:
        if s["status"] == "running":
            s["status"] = "done"
            s["elapsed"] = _format_elapsed(time.time() - s["start_time"])

    _steps.append({
        "name": description,
        "status": "running",
        "result": "",
        "elapsed": "",
        "start_time": time.time(),
    })
    _step_status = ""
    _refresh_live()


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
    """Context manager for spinner during subprocess operations.

    Updates the live display's transient status line.
    """
    global _step_status
    if _silent_mode:
        yield
        return

    if _live is not None:
        # Show as transient status under current step
        _step_status = description
        _refresh_live()
        yield
        _step_status = ""
        _refresh_live()
    else:
        # Standalone — print inline
        console.print(f"  [dim]\u23f3 {description}[/dim]", end="")
        yield
        console.print(" [done]\u2714[/done]")
