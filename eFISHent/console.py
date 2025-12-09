"""Rich console utilities for eFISHent CLI."""

from contextlib import contextmanager
from typing import Optional

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


def print_completion(duration: str) -> None:
    """Print the completion message."""
    if _silent_mode:
        return
    console.print(
        Panel(
            f"[success]eFISHent has finished successfully![/success]\n"
            f"[dim]Completed in {duration}[/dim]",
            border_style="green",
        )
    )


def print_stage(order: int, total: int, description: str) -> None:
    """Print a pipeline stage indicator with clear header."""
    global _current_stage, _total_stages
    if _silent_mode:
        return

    _current_stage = order
    _total_stages = total

    # Print a clear stage header
    console.print()
    console.print(f"[stage]Step {order}/{total}:[/stage] [current]{description}[/current]")

    # Update progress bar if active
    if _pipeline_progress is not None and _pipeline_task_id is not None:
        _pipeline_progress.update(
            _pipeline_task_id,
            completed=order - 1,  # Progress shows completed stages
            description=f"Running step {order}/{total}...",
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
            description=f"Completed step {_current_stage}/{_total_stages}",
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


def print_probe_table(df: pd.DataFrame, max_rows: int = 10) -> None:
    """Display probe data as a Rich table."""
    if _silent_mode:
        return

    table = Table(
        title="Probe Summary", show_header=True, header_style="bold cyan", box=None
    )
    table.add_column("Name", style="green")
    table.add_column("Sequence", style="dim", max_width=25)
    table.add_column("Length", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("TM", justify="right")
    table.add_column("deltaG", justify="right")

    for _, row in df.head(max_rows).iterrows():
        seq = row["sequence"]
        seq_display = seq[:22] + "..." if len(seq) > 25 else seq
        table.add_row(
            str(row["name"]),
            seq_display,
            str(row["length"]),
            f"{row['GC']:.1f}",
            f"{row['TM']:.1f}",
            f"{row['deltaG']:.1f}",
        )

    if len(df) > max_rows:
        table.add_row("...", f"({len(df) - max_rows} more)", "", "", "", "")

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
