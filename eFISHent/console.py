"""Rich console utilities for eFISHent CLI."""

import logging
from contextlib import contextmanager

from rich.panel import Panel
from rich.table import Table

from . import pipeline_display
from .completion_report import (
    _build_verification_summary as _build_verification_summary,
)
from .completion_report import print_completion as print_completion
from .console_base import THEME as THEME
from .console_base import console as console
from .console_base import is_silent as is_silent
from .console_base import set_silent_mode as set_silent_mode
from .filtering_funnel import _build_funnel_table as _build_funnel_table
from .filtering_funnel import _get_drop_pct as _get_drop_pct
from .filtering_funnel import _get_stage_color as _get_stage_color
from .filtering_funnel import get_funnel_data as get_funnel_data
from .filtering_funnel import print_filtering_funnel as print_filtering_funnel
from .filtering_funnel import record_funnel_stage as record_funnel_stage
from .filtering_funnel import reset_funnel_data as reset_funnel_data
from .pipeline_display import _refresh_live
from .pipeline_display import add_cached_step as add_cached_step
from .pipeline_display import pipeline_progress as pipeline_progress
from .pipeline_display import print_analysis_stage as print_analysis_stage
from .pipeline_display import print_candidate_count as print_candidate_count
from .pipeline_display import print_stage as print_stage
from .probe_display import _build_coverage_map as _build_coverage_map
from .probe_display import _build_probe_table as _build_probe_table
from .probe_display import _compute_probe_quality as _compute_probe_quality
from .probe_display import print_probe_table as print_probe_table


class LiveLogHandler(logging.Handler):
    """Logging handler that routes messages to the live step display."""

    def emit(self, record: logging.LogRecord) -> None:
        if pipeline_display._live is not None and not is_silent():
            # Show the message as transient status under current step
            msg = record.getMessage()
            # Strip leading whitespace for indented messages (like transcript list)
            msg = msg.strip()
            if msg:
                pipeline_display._step_status = msg
                _refresh_live()
        elif not is_silent():
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
    if is_silent():
        return
    from rich.rule import Rule

    console.print()
    console.print(f"  [bold cyan]eFISHent v{version}[/bold cyan]  [dim]RNA FISH probe designer[/dim]")
    console.print(Rule(style="blue"))


def print_parameter_warnings(warnings: list[str]) -> None:
    """Print parameter validation warnings in a panel."""
    if is_silent() or not warnings:
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
    if is_silent():
        return
    # If live display is active, show as transient status
    if pipeline_display._live is not None:
        pipeline_display._step_status = f"\u26a0 {message}"
        _refresh_live()
    else:
        console.print(f"[warning]Warning:[/warning] {message}")


def print_error_panel(title: str, message: str, hint: str | None = None) -> None:
    """Display an error in a formatted panel."""
    # Stop live display before showing error so it doesn't interfere
    if pipeline_display._live is not None:
        pipeline_display._live.stop()
    content = f"[error]{message}[/error]"
    if hint:
        content += f"\n\n[dim]Hint: {hint}[/dim]"
    console.print(Panel(content, title=f"[error]{title}[/error]", border_style="red"))


def print_dependency_check(results: dict[str, dict]) -> None:
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


def print_missing_deps_error(missing: list[str]) -> None:
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
    if is_silent():
        yield
        return

    if pipeline_display._live is not None:
        # Show as transient status under current step
        pipeline_display._step_status = description
        _refresh_live()
        yield
        pipeline_display._step_status = ""
        _refresh_live()
    else:
        # Standalone — print inline
        console.print(f"  [dim]\u23f3 {description}[/dim]", end="")
        yield
        console.print(" [done]\u2714[/done]")
