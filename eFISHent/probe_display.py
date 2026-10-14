"""Rich renderables for designed probes: probe table and gene coverage map."""

import pandas as pd
from rich.table import Table
from rich.text import Text

from .console_base import console, is_silent


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
    if is_silent():
        return
    table = _build_probe_table(df, max_rows)
    console.print(table)
    console.print()


def _build_coverage_map(
    probe_df: pd.DataFrame,
    gene_length: int,
    coverage_pct: float,
) -> Text | None:
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
