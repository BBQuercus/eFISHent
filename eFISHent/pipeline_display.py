"""Live pipeline step display: step tracking and its Rich rendering."""

import time
from contextlib import contextmanager

from rich.console import Console
from rich.live import Live
from rich.table import Table

from .console_base import _format_elapsed, console, is_silent
from .filtering_funnel import _funnel_data

# Module-level state for progress tracking
_current_stage: int = 0
_total_stages: int = 8

# Live display state
_live: Live | None = None
_steps: list[dict] = []  # Each: {name, status, result, elapsed, start_time}
_step_status: str = ""  # Transient status line for current step
_pipeline_mode: str = "pipeline"

# Braille spinner frames
_SPINNER_FRAMES = "\u280b\u2819\u2839\u2838\u283c\u2834\u2826\u2827\u2807\u280f"


def _step_row_style(step: dict, spinner_char: str) -> tuple[str, str, str, str]:
    """Return (icon, elapsed, name_style, result_style) for a pipeline step."""
    status = step["status"]
    if status == "done":
        return "[green]\u2714[/green]", step.get("elapsed", ""), "dim", "dim"
    if status == "running":
        # Compute elapsed live
        elapsed = _format_elapsed(time.time() - step["start_time"])
        return f"[cyan]{spinner_char}[/cyan]", elapsed, "bold", "cyan"
    return " ", "", "dim", "dim"


def _add_step_rows(table: Table, step: dict, spinner_char: str) -> None:
    """Add a pipeline step (and its transient status, if running) to the table."""
    name = step["name"]
    result = step.get("result", "")
    icon, elapsed, name_style, result_style = _step_row_style(step, spinner_char)

    table.add_row(
        icon,
        f"[{name_style}]{name}[/{name_style}]",
        "",
        f"[{result_style}]{result}[/{result_style}]" if result else "",
        f"[dim]{elapsed}[/dim]" if elapsed else "",
    )

    # Show transient status under the running step
    if step["status"] == "running" and _step_status:
        table.add_row(
            "",
            f"  [cyan]{spinner_char}[/cyan] [dim]{_step_status}[/dim]",
            "",
            "",
            "",
        )


def _add_total_elapsed_rows(table: Table) -> None:
    """Add the total elapsed time since the first step to the bottom of the table."""
    if not _steps:
        return
    first_start = _steps[0].get("start_time")
    if first_start is None:
        return
    total = _format_elapsed(time.time() - first_start)
    table.add_row("", "", "", "", "")
    table.add_row(
        "",
        "",
        "",
        f"[dim]{total} elapsed[/dim]",
        "",
    )


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
            has_running = has_running or step["status"] == "running"
            _add_step_rows(table, step, spinner_char)

        # Show status line when no step is running (e.g. dependency tasks)
        if not has_running and _step_status:
            table.add_row(
                f"[cyan]{spinner_char}[/cyan]",
                f"[dim]{_step_status}[/dim]",
                "",
                "",
                "",
            )

        _add_total_elapsed_rows(table)

        yield table


# Singleton display instance — reused across pipeline runs
_pipeline_display = _PipelineDisplay()


def _refresh_live() -> None:
    """Force a refresh of the live display if active."""
    if _live is not None:
        _live.refresh()


def add_cached_step(order: int, total: int, description: str) -> None:
    """Add a step that was already completed (cached by Luigi)."""
    if is_silent():
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
    if is_silent():
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
    if is_silent():
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

    if is_silent():
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
    if is_silent():
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
