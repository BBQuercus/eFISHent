"""Matplotlib panels for the probe set analysis PDF."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable


def histplot(ax: plt.Axes, data: list[int], title: str, min_value: int):
    """Basic narrow count based histogram."""
    ax.set_title(title)
    if not data or all(v == data[0] for v in data):
        ax.hist(data, bins=10)
    else:
        bins = np.arange(min_value, max(data) + 4) + 0.5
        ax.hist(data, bins=bins)
        ax.set(xticks=bins[1:] - 0.5)
    ax.set(ylabel="Count")


def boxplot(ax: plt.Axes, data: list[int | float], title: str, ylabel: str):
    """Boxplot with kde-esque lines next to it."""
    ax.set_title(title)
    ax.boxplot(data)
    ax.eventplot(
        data,
        orientation="vertical",
        linewidths=0.2,
        lineoffsets=0.8,
        linelengths=0.1,
    )
    ax.set(xticks=[], ylabel=ylabel)


def matrix(ax: plt.Axes, data: np.ndarray, title: str, fig):
    """Image with colorbar in same axis."""
    ax.set_title(title)
    image = ax.imshow(data, vmin=0, vmax=1, cmap="coolwarm")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(image, cax=cax, orientation="vertical")


def plot_probe_coverage(ax: plt.Axes, df: pd.DataFrame, gene_len: int):
    """Gene coverage track with probe positions and large gaps annotated."""
    coverage = np.zeros(gene_len)
    mapped = df[df["start"] >= 0]

    for _, row in mapped.iterrows():
        s, e = int(row["start"]), int(row["end"])
        coverage[s:e] = 1

    ax.set_title("Gene coverage")
    ax.fill_between(range(gene_len), coverage, alpha=0.6, color="steelblue")

    # Mark probe positions
    for _, row in mapped.iterrows():
        mid = (row["start"] + row["end"]) / 2
        ax.axvline(mid, color="steelblue", alpha=0.3, linewidth=0.5)

    ax.set(ylabel="Coverage", xlabel="Position (nt)", ylim=(-0.1, 1.5))
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["No", "Yes"])

    # Add gap annotations
    if len(mapped) > 1:
        sorted_mapped = mapped.sort_values("start")
        prev_end = 0
        gaps = []
        for _, row in sorted_mapped.iterrows():
            gap = row["start"] - prev_end
            if gap > 100:
                gaps.append((prev_end, int(row["start"]), gap))
            prev_end = max(prev_end, int(row["end"]))
        for gs, ge, size in gaps[:5]:  # Show max 5 largest gaps
            mid = (gs + ge) / 2
            ax.annotate(
                f"{size}nt gap",
                xy=(mid, 1.1), fontsize=6, ha="center", color="red", alpha=0.7,
            )
