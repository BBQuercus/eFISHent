"""Tests for console UX features."""

import io

import pandas as pd
import pytest
from rich.console import Console

from eFISHent.console import (
    get_funnel_data,
    print_completion,
    print_filtering_funnel,
    print_probe_table,
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


class TestFunnelData:
    """Tests for filtering funnel data accumulation."""

    def test_record_and_get(self):
        record_funnel_stage("Generated", 1000)
        record_funnel_stage("TM/GC filtered", 800)
        data = get_funnel_data()
        assert len(data) == 2
        assert data[0] == ("Generated", 1000)
        assert data[1] == ("TM/GC filtered", 800)

    def test_reset(self):
        record_funnel_stage("Generated", 1000)
        reset_funnel_data()
        assert get_funnel_data() == []

    def test_get_returns_copy(self):
        record_funnel_stage("Generated", 1000)
        data = get_funnel_data()
        data.append(("fake", 0))
        assert len(get_funnel_data()) == 1


class TestPrintFilteringFunnel:
    """Tests for the filtering funnel visualization."""

    def test_funnel_with_data(self, capsys):
        record_funnel_stage("Generated", 1247)
        record_funnel_stage("TM/GC content", 892)
        record_funnel_stage("Genome alignment", 583)
        record_funnel_stage("K-mer frequency", 541)
        record_funnel_stage("Secondary structure", 498)
        record_funnel_stage("Optimization", 48)
        print_filtering_funnel()
        captured = capsys.readouterr()
        assert "Filtering Funnel" in captured.out
        assert "1,247" in captured.out
        assert "48" in captured.out

    def test_funnel_empty(self, capsys):
        print_filtering_funnel()
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_funnel_silent_mode(self, capsys):
        set_silent_mode(True)
        record_funnel_stage("Generated", 1000)
        print_filtering_funnel()
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_funnel_single_stage(self, capsys):
        record_funnel_stage("Generated", 500)
        print_filtering_funnel()
        captured = capsys.readouterr()
        assert "500" in captured.out


class TestPrintCompletion:
    """Tests for the completion summary panel."""

    def test_completion_with_summary(self, capsys):
        summary = {
            "gene_name": "aad4",
            "probe_count": 48,
            "initial_count": 1247,
            "coverage_pct": 78.3,
            "length_range": (20, 25),
            "length_median": 22,
            "tm_range": (41.2, 58.7),
            "tm_median": 49.3,
            "gc_range": (24.0, 68.2),
            "gc_median": 45.1,
        }
        print_completion("2m 15s", ["/tmp/out.fasta"], summary)
        captured = capsys.readouterr()
        assert "aad4" in captured.out
        assert "48" in captured.out
        assert "1,247" in captured.out
        assert "78.3%" in captured.out
        assert "20-25" in captured.out
        assert "/tmp/out.fasta" in captured.out
        assert "2m 15s" in captured.out

    def test_completion_without_summary(self, capsys):
        print_completion("1.5s", ["/tmp/out.fasta"])
        captured = capsys.readouterr()
        assert "Design Complete" in captured.out
        assert "/tmp/out.fasta" in captured.out

    def test_completion_silent_mode(self, capsys):
        set_silent_mode(True)
        print_completion("1s", ["/tmp/out.fasta"])
        captured = capsys.readouterr()
        assert captured.out == ""


class TestPrintProbeTable:
    """Tests for the probe table display."""

    @pytest.fixture
    def sample_df(self):
        return pd.DataFrame(
            {
                "name": [f"probe-{i}" for i in range(1, 21)],
                "sequence": ["ATCGATCGATCGATCGATCG"] * 20,
                "length": [20] * 20,
                "GC": [50.0] * 20,
                "TM": [55.0] * 20,
                "deltaG": [-5.0] * 20,
                "kmers": [3] * 20,
                "count": [0] * 20,
            }
        )

    def test_shows_all_probes_default(self, sample_df, capsys):
        """Default max_rows=100 should show all 20 probes."""
        print_probe_table(sample_df)
        captured = capsys.readouterr()
        assert "probe-1" in captured.out
        assert "probe-20" in captured.out
        assert "more" not in captured.out

    def test_silent_mode(self, sample_df, capsys):
        set_silent_mode(True)
        print_probe_table(sample_df)
        captured = capsys.readouterr()
        assert captured.out == ""
