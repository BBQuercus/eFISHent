"""Tests for quality scores, recommendation, and off-target cap."""

import re
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from eFISHent.cleanup import CleanUpOutput


@pytest.fixture
def perfect_probe_df():
    """DataFrame with a near-perfect probe."""
    return pd.DataFrame(
        {
            "name": ["probe-1"],
            "sequence": ["ATGCATGCATGCATGCATGCATGC"],
            "length": [24],
            "start": [0],
            "end": [24],
            "GC": [50.0],
            "TM": [55.0],
            "deltaG": [0.0],
            "kmers": [0],
            "count": [0],
            "txome_off_targets": [0],
            "off_target_genes": [""],
            "worst_match": [""],
            "expression_risk": [""],
        }
    )


@pytest.fixture
def bad_probe_df():
    """DataFrame with a poor-quality probe."""
    return pd.DataFrame(
        {
            "name": ["probe-bad"],
            "sequence": ["GCGCGCGCGCGCGCGCGCGCGCGC"],
            "length": [24],
            "start": [0],
            "end": [24],
            "GC": [100.0],
            "TM": [80.0],
            "deltaG": [-15.0],
            "kmers": [10],
            "count": [5],
            "txome_off_targets": [10],
            "off_target_genes": ["TP53(5), MYC(5)"],
            "worst_match": ["100%/24bp/0mm"],
            "expression_risk": ["TP53:HIGH(200)"],
        }
    )


@pytest.fixture
def multi_probe_df():
    """DataFrame with multiple probes of varying quality."""
    return pd.DataFrame(
        {
            "name": ["probe-1", "probe-2", "probe-3", "probe-4"],
            "sequence": [
                "ATGCATGCATGCATGCATGCATGC",
                "GCGCGCGCGCGCGCGCGCGCGCGC",
                "ATATATATATATATATATATATATAT",
                "ATGCATGCATGCATGCATGCATGC",
            ],
            "length": [24, 24, 24, 24],
            "start": [0, 24, 48, 72],
            "end": [24, 48, 72, 96],
            "GC": [50.0, 100.0, 0.0, 50.0],
            "TM": [55.0, 80.0, 30.0, 55.0],
            "deltaG": [0.0, -15.0, -12.0, -2.0],
            "kmers": [0, 10, 8, 1],
            "count": [0, 5, 3, 0],
            "txome_off_targets": [0, 10, 5, 0],
            "off_target_genes": [
                "",
                "TP53(5), MYC(5)",
                "ACTB(3), GAPDH(2)",
                "",
            ],
            "worst_match": ["", "100%/24bp/0mm", "95%/20bp/1mm", ""],
            "expression_risk": ["", "TP53:HIGH(200)", "", ""],
        }
    )


class TestComputeQualityScores:
    def test_perfect_probe_high_score(self, perfect_probe_df):
        """A perfect probe should score near 100."""
        scores = CleanUpOutput._compute_quality_scores(perfect_probe_df)
        assert scores.iloc[0] >= 90.0

    def test_bad_probe_low_score(self, bad_probe_df):
        """A bad probe should score low."""
        scores = CleanUpOutput._compute_quality_scores(bad_probe_df)
        assert scores.iloc[0] < 50.0

    def test_scores_in_range(self, multi_probe_df):
        """All scores should be in [0, 100]."""
        scores = CleanUpOutput._compute_quality_scores(multi_probe_df)
        assert (scores >= 0).all()
        assert (scores <= 100).all()

    def test_good_probe_beats_bad_probe(self, multi_probe_df):
        """Probe with good stats should score higher than probe with bad stats."""
        scores = CleanUpOutput._compute_quality_scores(multi_probe_df)
        # probe-1 (idx 0) has ideal stats; probe-2 (idx 1) has terrible stats
        assert scores.iloc[0] > scores.iloc[1]

    def test_score_is_series(self, multi_probe_df):
        """Output should be a pandas Series with correct length."""
        scores = CleanUpOutput._compute_quality_scores(multi_probe_df)
        assert isinstance(scores, pd.Series)
        assert len(scores) == len(multi_probe_df)

    def test_zero_off_targets_boosts_score(self):
        """Probes with zero off-targets should score better than those with many."""
        df = pd.DataFrame(
            {
                "GC": [50.0, 50.0],
                "TM": [55.0, 55.0],
                "deltaG": [0.0, 0.0],
                "kmers": [0, 0],
                "count": [0, 5],
                "txome_off_targets": [0, 5],
            }
        )
        scores = CleanUpOutput._compute_quality_scores(df)
        assert scores.iloc[0] > scores.iloc[1]


class TestComputeRecommendation:
    def test_pass_high_quality_no_off_targets(self):
        """High quality with no off-targets should be PASS."""
        row = pd.Series(
            {"quality": 85.0, "txome_off_targets": 0, "expression_risk": ""}
        )
        assert CleanUpOutput._compute_recommendation(row) == "PASS"

    def test_flag_low_quality(self):
        """Quality < 70 should be FLAG."""
        row = pd.Series(
            {"quality": 55.0, "txome_off_targets": 0, "expression_risk": ""}
        )
        result = CleanUpOutput._compute_recommendation(row)
        assert result.startswith("FLAG")
        assert "low_quality" in result

    def test_flag_many_off_targets(self):
        """Many transcriptome off-targets (>3) should be FLAG."""
        row = pd.Series(
            {"quality": 80.0, "txome_off_targets": 5, "expression_risk": ""}
        )
        result = CleanUpOutput._compute_recommendation(row)
        assert result.startswith("FLAG")
        assert "many_off_targets" in result

    def test_flag_high_expression_risk(self):
        """HIGH expression risk should be FLAG."""
        row = pd.Series(
            {
                "quality": 80.0,
                "txome_off_targets": 0,
                "expression_risk": "TP53:HIGH(200)",
            }
        )
        result = CleanUpOutput._compute_recommendation(row)
        assert result.startswith("FLAG")
        assert "high_expression_risk" in result

    def test_fail_very_low_quality(self):
        """Quality < 30 should be FAIL."""
        row = pd.Series(
            {"quality": 20.0, "txome_off_targets": 0, "expression_risk": ""}
        )
        assert CleanUpOutput._compute_recommendation(row) == "FAIL"

    def test_fail_takes_priority_over_flag(self):
        """FAIL (quality<30) should take priority even if FLAG conditions exist."""
        row = pd.Series(
            {
                "quality": 15.0,
                "txome_off_targets": 10,
                "expression_risk": "TP53:HIGH(500)",
            }
        )
        assert CleanUpOutput._compute_recommendation(row) == "FAIL"

    def test_multiple_flags(self):
        """Multiple FLAG conditions should all appear."""
        row = pd.Series(
            {
                "quality": 50.0,
                "txome_off_targets": 5,
                "expression_risk": "TP53:HIGH(200)",
            }
        )
        result = CleanUpOutput._compute_recommendation(row)
        assert result.startswith("FLAG")
        assert "low_quality" in result
        assert "many_off_targets" in result
        assert "high_expression_risk" in result

    def test_pass_boundary(self):
        """Quality exactly 70 with no issues should be PASS."""
        row = pd.Series(
            {"quality": 70.0, "txome_off_targets": 0, "expression_risk": ""}
        )
        assert CleanUpOutput._compute_recommendation(row) == "PASS"

    def test_fail_boundary(self):
        """Quality exactly 30 should NOT be FAIL (< 30, not <=)."""
        row = pd.Series(
            {"quality": 30.0, "txome_off_targets": 0, "expression_risk": ""}
        )
        result = CleanUpOutput._compute_recommendation(row)
        assert result != "FAIL"


class TestApplyOffTargetCap:
    def _make_df_with_off_targets(self, off_target_genes, quality_scores):
        """Helper to build a DataFrame for off-target cap testing."""
        n = len(off_target_genes)
        return pd.DataFrame(
            {
                "name": [f"probe-{i+1}" for i in range(n)],
                "sequence": ["ATGCATGCATGCATGCATGCATGC"] * n,
                "length": [24] * n,
                "start": list(range(0, n * 24, 24)),
                "end": list(range(24, (n + 1) * 24, 24)),
                "GC": [50.0] * n,
                "TM": [55.0] * n,
                "deltaG": [0.0] * n,
                "kmers": [0] * n,
                "count": [0] * n,
                "txome_off_targets": [1 if g else 0 for g in off_target_genes],
                "off_target_genes": off_target_genes,
                "worst_match": ["" if not g else "90%/20bp/2mm" for g in off_target_genes],
                "expression_risk": [""] * n,
                "quality": quality_scores,
            }
        )

    def _mock_probe_config(self, cap_value):
        """Create a mock ProbeConfig with the given cap value."""
        mock_config = MagicMock()
        mock_config.max_probes_per_off_target = cap_value
        return mock_config

    def test_cap_disabled(self):
        """When cap=0, no probes should be removed."""
        task = CleanUpOutput()
        df = self._make_df_with_off_targets(
            ["TP53(1)", "TP53(1)", "TP53(1)"],
            [80.0, 60.0, 70.0],
        )
        with patch("eFISHent.cleanup.ProbeConfig", return_value=self._mock_probe_config(0)):
            result = task._apply_off_target_cap(df)
            assert len(result) == 3

    def test_cap_removes_lowest_quality(self):
        """When cap is exceeded, the lowest-quality probe should be removed."""
        task = CleanUpOutput()
        # 3 probes all hitting TP53, cap=2
        df = self._make_df_with_off_targets(
            ["TP53(1)", "TP53(1)", "TP53(1)"],
            [80.0, 50.0, 70.0],
        )
        with patch("eFISHent.cleanup.ProbeConfig", return_value=self._mock_probe_config(2)):
            result = task._apply_off_target_cap(df)
            assert len(result) == 2
            # probe-2 (quality=50) should have been removed
            remaining_qualities = sorted(result["quality"].tolist())
            assert 50.0 not in remaining_qualities

    def test_cap_no_removal_when_under(self):
        """When under cap, no probes should be removed."""
        task = CleanUpOutput()
        # 2 probes hitting TP53, cap=3
        df = self._make_df_with_off_targets(
            ["TP53(1)", "TP53(1)", ""],
            [80.0, 70.0, 90.0],
        )
        with patch("eFISHent.cleanup.ProbeConfig", return_value=self._mock_probe_config(3)):
            result = task._apply_off_target_cap(df)
            assert len(result) == 3

    def test_cap_no_off_target_column(self):
        """When off_target_genes column is missing, should return unchanged."""
        task = CleanUpOutput()
        df = pd.DataFrame(
            {
                "name": ["probe-1"],
                "quality": [80.0],
            }
        )
        with patch("eFISHent.cleanup.ProbeConfig", return_value=self._mock_probe_config(2)):
            result = task._apply_off_target_cap(df)
            assert len(result) == 1
