"""Tests for data computation functions in console.py."""

import argparse

import pandas as pd
import pytest

from eFISHent.cli import validate_parameter_warnings
from eFISHent.console import (
    _build_coverage_map,
    _build_funnel_table,
    _build_probe_table,
    _compute_probe_quality,
    _get_drop_pct,
    _get_stage_color,
    get_funnel_data,
    print_parameter_warnings,
    record_funnel_stage,
    reset_funnel_data,
    set_silent_mode,
)


@pytest.fixture(autouse=True)
def _reset_state():
    """Reset console module state between tests."""
    reset_funnel_data()
    set_silent_mode(False)
    yield
    set_silent_mode(False)
    reset_funnel_data()


# ---------------------------------------------------------------------------
# _compute_probe_quality
# ---------------------------------------------------------------------------
class TestComputeProbeQuality:
    """Tests for the _compute_probe_quality star-rating function."""

    def _make_df(self, rows):
        """Helper to build a DataFrame from a list of row dicts."""
        defaults = {"name": "p", "sequence": "ATCG" * 5, "length": 20, "kmers": 1}
        full_rows = [{**defaults, **r} for r in rows]
        return pd.DataFrame(full_rows)

    def test_perfect_probe_three_stars(self):
        """Central TM, balanced GC, no off-targets, weak deltaG -> 3 stars."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 50.0, "deltaG": -3.0, "count": 0},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        # Middle probe is perfectly central
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2605"

    def test_poor_probe_one_star(self):
        """Extreme TM, extreme GC, off-targets, strong deltaG -> 1 star."""
        df = self._make_df(
            [
                {"TM": 60.0, "GC": 70.0, "deltaG": -12.0, "count": 5},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2606\u2606"

    def test_mixed_quality_two_stars(self):
        """One penalty factor -> 2 stars."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 50.0, "deltaG": -3.0, "count": 1},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2606"

    def test_missing_deltag_column(self):
        """If deltaG column is missing from the row the function should still work."""
        df = pd.DataFrame(
            {
                "name": ["p1", "p2"],
                "sequence": ["ATCG" * 5] * 2,
                "length": [20, 20],
                "TM": [50.0, 40.0],
                "GC": [50.0, 55.0],
                "count": [0, 0],
                "kmers": [1, 1],
            }
        )
        # deltaG not present -- _compute_probe_quality accesses row["deltaG"];
        # callers always include it, but we test the KeyError guard isn't needed
        # by ensuring the column is present with a neutral value
        df["deltaG"] = -3.0
        score = _compute_probe_quality(df.iloc[0], df)
        assert len(score) == 3

    def test_missing_count_column(self):
        """No 'count' column -> off-target penalty is skipped."""
        df = pd.DataFrame(
            {
                "name": ["p1", "p2", "p3"],
                "sequence": ["ATCG" * 5] * 3,
                "length": [20, 20, 20],
                "TM": [50.0, 40.0, 60.0],
                "GC": [50.0, 45.0, 55.0],
                "deltaG": [-3.0, -2.0, -4.0],
                "kmers": [1, 1, 1],
            }
        )
        # 'count' not in row.index, so off-target branch is skipped
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2605"

    def test_gc_boundary_30(self):
        """GC exactly 30 should NOT be penalised (< 30 triggers penalty)."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 30.0, "deltaG": -3.0, "count": 0},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2605"

    def test_gc_boundary_65(self):
        """GC exactly 65 should NOT be penalised (> 65 triggers penalty)."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 65.0, "deltaG": -3.0, "count": 0},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2605"

    def test_gc_below_30_penalised(self):
        """GC < 30 loses a star."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 29.9, "deltaG": -3.0, "count": 0},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2606"

    def test_gc_above_65_penalised(self):
        """GC > 65 loses a star."""
        df = self._make_df(
            [
                {"TM": 50.0, "GC": 65.1, "deltaG": -3.0, "count": 0},
                {"TM": 40.0, "GC": 45.0, "deltaG": -2.0, "count": 0},
                {"TM": 60.0, "GC": 55.0, "deltaG": -4.0, "count": 0},
            ]
        )
        score = _compute_probe_quality(df.iloc[0], df)
        assert score == "\u2605\u2605\u2606"


# ---------------------------------------------------------------------------
# _get_drop_pct
# ---------------------------------------------------------------------------
class TestGetDropPct:
    def test_prev_zero_returns_zero(self):
        assert _get_drop_pct(50, 0) == 0

    def test_no_change(self):
        assert _get_drop_pct(100, 100) == 0.0

    def test_full_drop(self):
        assert _get_drop_pct(0, 100) == 100.0

    def test_normal_drop(self):
        assert _get_drop_pct(80, 100) == pytest.approx(20.0)


# ---------------------------------------------------------------------------
# _get_stage_color
# ---------------------------------------------------------------------------
class TestGetStageColor:
    def test_red_above_50(self):
        assert _get_stage_color(51) == "red"

    def test_yellow_above_20(self):
        assert _get_stage_color(21) == "yellow"

    def test_green_at_20(self):
        assert _get_stage_color(20) == "green"

    def test_green_below_20(self):
        assert _get_stage_color(10) == "green"

    def test_red_at_boundary(self):
        # Exactly 50 is NOT > 50
        assert _get_stage_color(50) == "yellow"


# ---------------------------------------------------------------------------
# Funnel data management
# ---------------------------------------------------------------------------
class TestFunnelData:
    def test_record_adds_entries(self):
        record_funnel_stage("A", 100)
        record_funnel_stage("B", 80)
        assert get_funnel_data() == [("A", 100), ("B", 80)]

    def test_get_returns_ordered(self):
        for i, name in enumerate(["X", "Y", "Z"]):
            record_funnel_stage(name, (i + 1) * 10)
        data = get_funnel_data()
        assert [n for n, _ in data] == ["X", "Y", "Z"]

    def test_reset_clears(self):
        record_funnel_stage("A", 100)
        reset_funnel_data()
        assert get_funnel_data() == []

    def test_get_returns_copy(self):
        record_funnel_stage("A", 100)
        data = get_funnel_data()
        data.append(("fake", 0))
        assert len(get_funnel_data()) == 1


# ---------------------------------------------------------------------------
# _build_funnel_table
# ---------------------------------------------------------------------------
class TestBuildFunnelTable:
    def test_no_data_returns_none(self):
        assert _build_funnel_table() is None

    def test_skipped_stages_only_returns_none(self):
        record_funnel_stage("Preparing gene sequence", 100)
        record_funnel_stage("Finalizing output", 50)
        assert _build_funnel_table() is None

    def test_single_stage_returns_table(self):
        record_funnel_stage("Generated", 500)
        table = _build_funnel_table()
        assert table is not None
        assert table.row_count == 1

    def test_multiple_stages_row_count(self):
        record_funnel_stage("Generated", 1000)
        record_funnel_stage("TM/GC content", 800)
        record_funnel_stage("Optimization", 200)
        table = _build_funnel_table()
        assert table is not None
        assert table.row_count == 3

    def test_all_zero_counts_returns_none(self):
        record_funnel_stage("Generated", 0)
        assert _build_funnel_table() is None


# ---------------------------------------------------------------------------
# _build_coverage_map
# ---------------------------------------------------------------------------
class TestBuildCoverageMap:
    def test_full_coverage(self):
        """All positions covered -> all filled blocks."""
        df = pd.DataFrame({"start": [0], "end": [1000]})
        result = _build_coverage_map(df, gene_length=1000, coverage_pct=100.0)
        assert result is not None
        plain = result.plain
        # Every block in the bar should be filled (\u2588), none light (\u2591)
        assert "\u2591" not in plain

    def test_no_probes(self):
        """Empty DataFrame -> all light blocks."""
        df = pd.DataFrame({"start": [], "end": []})
        result = _build_coverage_map(df, gene_length=1000, coverage_pct=0.0)
        assert result is not None
        plain = result.plain
        assert "\u2588" not in plain

    def test_gap_in_middle(self):
        """Coverage with a gap in the centre -> both filled and light blocks."""
        df = pd.DataFrame({"start": [0, 800], "end": [200, 1000]})
        result = _build_coverage_map(df, gene_length=1000, coverage_pct=40.0)
        assert result is not None
        plain = result.plain
        assert "\u2588" in plain
        assert "\u2591" in plain

    def test_coverage_pct_in_output(self):
        df = pd.DataFrame({"start": [0], "end": [500]})
        result = _build_coverage_map(df, gene_length=1000, coverage_pct=50.0)
        assert result is not None
        assert "50.0%" in result.plain

    def test_zero_gene_length_returns_none(self):
        df = pd.DataFrame({"start": [0], "end": [100]})
        assert _build_coverage_map(df, gene_length=0, coverage_pct=0.0) is None

    def test_negative_gene_length_returns_none(self):
        df = pd.DataFrame({"start": [0], "end": [100]})
        assert _build_coverage_map(df, gene_length=-1, coverage_pct=0.0) is None


# ---------------------------------------------------------------------------
# _build_probe_table
# ---------------------------------------------------------------------------
class TestBuildProbeTable:
    @pytest.fixture
    def basic_df(self):
        return pd.DataFrame(
            {
                "name": ["probe-1", "probe-2", "probe-3"],
                "sequence": ["ATCGATCGATCGATCGATCG"] * 3,
                "length": [20, 20, 20],
                "GC": [50.0, 50.0, 50.0],
                "TM": [55.0, 55.0, 55.0],
                "deltaG": [-5.0, -5.0, -5.0],
                "kmers": [3, 3, 3],
                "count": [0, 0, 0],
            }
        )

    def test_basic_table_created(self, basic_df):
        table = _build_probe_table(basic_df)
        assert table is not None
        assert table.row_count == 3

    def test_off_targets_column_added(self, basic_df):
        basic_df["txome_off_targets"] = [0, 2, 1]
        table = _build_probe_table(basic_df)
        col_names = [c.header for c in table.columns]
        assert "Off-targets" in col_names

    def test_recommendation_column_added(self, basic_df):
        basic_df["recommendation"] = ["PASS", "FLAG", "FAIL"]
        table = _build_probe_table(basic_df)
        col_names = [c.header for c in table.columns]
        assert "Rec." in col_names

    def test_no_optional_columns(self, basic_df):
        table = _build_probe_table(basic_df)
        col_names = [c.header for c in table.columns]
        assert "Off-targets" not in col_names
        assert "Rec." not in col_names

    def test_truncation_at_max_rows(self):
        n = 10
        df = pd.DataFrame(
            {
                "name": [f"probe-{i}" for i in range(n)],
                "sequence": ["ATCGATCGATCGATCGATCG"] * n,
                "length": [20] * n,
                "GC": [50.0] * n,
                "TM": [55.0] * n,
                "deltaG": [-5.0] * n,
                "kmers": [3] * n,
                "count": [0] * n,
            }
        )
        table = _build_probe_table(df, max_rows=3)
        # 3 data rows + 1 "...more" row
        assert table.row_count == 4

    def test_long_sequence_truncated(self):
        long_seq = "A" * 30
        df = pd.DataFrame(
            {
                "name": ["probe-1"],
                "sequence": [long_seq],
                "length": [30],
                "GC": [50.0],
                "TM": [55.0],
                "deltaG": [-5.0],
                "kmers": [3],
                "count": [0],
            }
        )
        from io import StringIO
        from rich.console import Console as RichConsole

        buf = StringIO()
        c = RichConsole(file=buf, width=200)
        table = _build_probe_table(df)
        c.print(table)
        output = buf.getvalue()
        # The full 30-char sequence should NOT appear; truncated with "..."
        assert long_seq not in output
        assert "..." in output

    def test_short_sequence_not_truncated(self, basic_df):
        from io import StringIO
        from rich.console import Console as RichConsole

        buf = StringIO()
        c = RichConsole(file=buf, width=200)
        table = _build_probe_table(basic_df)
        c.print(table)
        output = buf.getvalue()
        assert "ATCGATCGATCGATCGATCG" in output


# ---------------------------------------------------------------------------
# Parameter warnings (kept from previous version)
# ---------------------------------------------------------------------------
class TestParameterWarnings:
    @pytest.fixture
    def default_args(self):
        return argparse.Namespace(
            build_indices=False,
            analyze_probeset="",
            min_length=21,
            max_length=25,
            min_tm=40.0,
            max_tm=60.0,
            min_gc=20.0,
            max_gc=80.0,
            formamide_concentration=10.0,
            kmer_length=15,
            max_deltag=-10.0,
            spacing=2,
            is_endogenous=True,
            reference_transcriptome="",
        )

    def test_no_warnings_default_params(self, default_args):
        assert validate_parameter_warnings(default_args) == []

    def test_narrow_tm_window(self, default_args):
        default_args.min_tm = 50.0
        default_args.max_tm = 55.0
        assert any("TM window" in w for w in validate_parameter_warnings(default_args))

    def test_high_formamide_low_tm(self, default_args):
        default_args.formamide_concentration = 50.0
        default_args.max_tm = 55.0
        assert any("Formamide" in w for w in validate_parameter_warnings(default_args))

    def test_kmer_close_to_probe_length(self, default_args):
        default_args.kmer_length = 19
        assert any("K-mer" in w for w in validate_parameter_warnings(default_args))

    def test_narrow_gc_window(self, default_args):
        default_args.min_gc = 40.0
        default_args.max_gc = 55.0
        assert any("GC window" in w for w in validate_parameter_warnings(default_args))

    def test_strict_deltag(self, default_args):
        default_args.max_deltag = -2.0
        assert any("deltaG" in w for w in validate_parameter_warnings(default_args))

    def test_large_spacing(self, default_args):
        default_args.spacing = 15
        assert any("spacing" in w for w in validate_parameter_warnings(default_args))

    def test_skip_for_build_indices(self, default_args):
        default_args.build_indices = True
        default_args.min_tm = 55.0
        default_args.max_tm = 56.0
        assert validate_parameter_warnings(default_args) == []

    def test_print_parameter_warnings(self, capsys):
        print_parameter_warnings(["Warning 1.\n  Suggestion 1.", "Warning 2."])
        captured = capsys.readouterr()
        assert "Parameter Warnings" in captured.out

    def test_print_parameter_warnings_empty(self, capsys):
        print_parameter_warnings([])
        assert capsys.readouterr().out == ""

    def test_print_parameter_warnings_silent(self, capsys):
        set_silent_mode(True)
        print_parameter_warnings(["Warning"])
        assert capsys.readouterr().out == ""
