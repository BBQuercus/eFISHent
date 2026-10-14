"""Swap Tm outliers in the selected probe set for better-matched alternatives."""

from collections.abc import Callable
from dataclasses import dataclass

import Bio.SeqRecord
import numpy as np
import pandas as pd

from .config import ProbeConfig
from .selection import is_overlapping


class TmUniformityMixin:
    """Refine an assigned probe set towards a uniform melting temperature."""

    def refine_tm_uniformity(
        self,
        assigned: list[str],
        sequences: list[Bio.SeqRecord.SeqRecord],
    ) -> list[str]:
        """Refine probe selection for Tm uniformity.

        For each probe whose Tm deviates > 5°C from the set median, check
        if an alternative probe at a nearby position has better Tm and swap.
        """
        from .basic_filtering import get_melting_temp

        if len(assigned) < 3:
            return assigned

        config = ProbeConfig()
        na = config.na_concentration
        formamide = config.formamide_concentration

        # Compute Tm for all assigned probes
        seq_map = {s.id: s for s in sequences}
        assigned_tms = {
            name: get_melting_temp(seq_map[name].seq, na, formamide)
            for name in assigned
            if name in seq_map
        }
        if not assigned_tms:
            return assigned

        assigned_set = set(assigned)
        search = _TmSwapSearch(
            unassigned=[s for s in sequences if s.id not in assigned_set],
            assigned=list(assigned),
            df=self.df,
            median_tm=float(np.median(list(assigned_tms.values()))),
            na=na,
            formamide=formamide,
            melting_temp=get_melting_temp,
        )

        swaps_made = 0
        for i, name in enumerate(search.assigned):
            if name not in assigned_tms:
                continue
            best_alt = self._tm_swap(search, i, assigned_tms[name])
            if best_alt is not None:
                search.assigned[i] = best_alt.id
                swaps_made += 1

        if swaps_made > 0:
            self.logger.debug(
                f"Tm uniformity: swapped {swaps_made} probes for better Tm match "
                f"(median Tm: {search.median_tm:.1f}°C)"
            )

        return search.assigned

    def _tm_swap(self, search: "_TmSwapSearch", index: int, tm: float):
        """Find a replacement for the assigned probe at index if its Tm is an outlier.

        Returns the alternative probe if it lowers the Tm deviation from the
        median by more than 1°C, otherwise None.
        """
        max_deviation = 5.0
        deviation = abs(tm - search.median_tm)
        if deviation <= max_deviation:
            return None

        name = search.assigned[index]
        probe_row = self.df[self.df["name"] == name]
        if probe_row.empty:
            return None

        probe_start = int(name.split("-")[-1])
        probe_end = int(probe_row.iloc[0]["end"])
        best_alt, best_alt_deviation = _find_best_tm_alternative(
            search, index, probe_start, probe_end, deviation
        )
        if best_alt and best_alt_deviation < deviation - 1.0:
            return best_alt
        return None


@dataclass
class _TmSwapSearch:
    """State shared while swapping Tm outliers for unassigned alternatives."""

    unassigned: list[Bio.SeqRecord.SeqRecord]
    assigned: list[str]
    df: pd.DataFrame
    median_tm: float
    na: float
    formamide: float
    melting_temp: Callable

    def deviation(self, seq) -> float:
        """Absolute difference of a sequence's Tm from the median Tm."""
        return abs(self.melting_temp(seq, self.na, self.formamide) - self.median_tm)


def _overlaps_other_assigned(
    alt_start: int,
    alt_end: int,
    new_assigned: list[str],
    skip_idx: int,
    df: pd.DataFrame,
) -> bool:
    """Check if an alternative probe overlaps any other assigned probe."""
    for j, other_name in enumerate(new_assigned):
        if j == skip_idx:
            continue
        other_row = df[df["name"] == other_name]
        if not other_row.empty:
            other_start = int(other_row.iloc[0]["start"])
            other_end = int(other_row.iloc[0]["end"])
            if is_overlapping((alt_start, alt_end), (other_start, other_end)):
                return True
    return False


def _find_best_tm_alternative(
    search: _TmSwapSearch,
    skip_idx: int,
    probe_start: int,
    probe_end: int,
    current_deviation: float,
):
    """Find the best alternative probe with lower Tm deviation."""
    best_alt = None
    best_alt_deviation = current_deviation

    for alt in search.unassigned:
        alt_start = int(alt.id.split("-")[-1])
        alt_end = alt_start + len(alt)

        # Alternative must be at a nearby position (overlapping region)
        if alt_end < probe_start or alt_start > probe_end:
            continue

        if _overlaps_other_assigned(
            alt_start, alt_end, search.assigned, skip_idx, search.df
        ):
            continue

        alt_deviation = search.deviation(alt.seq)
        if alt_deviation < best_alt_deviation:
            best_alt = alt
            best_alt_deviation = alt_deviation

    return best_alt, best_alt_deviation
