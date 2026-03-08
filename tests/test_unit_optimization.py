import shutil

import pandas as pd
import pytest

from eFISHent.optimization import OptimizeProbeCoverage
from eFISHent.optimization import fill_coverage_gaps
from eFISHent.optimization import greedy_model
from eFISHent.optimization import is_binding
from eFISHent.optimization import is_overlapping
from eFISHent.optimization import optimal_model

GLPK_AVAILABLE = shutil.which("glpsol") is not None


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"],
            "sequence": [
                "GTAATTACAAAATAAGCAACG",
                "GCTTGCTTTGAGATTTTGTTC",
                "TCAATTCTCTACTGTCTCAGT",
                "GTAATTACAAAATAAGCAACG",
                "GGGACGAATTCTTTGTCTATTC",
                "TGGATTTCATAATGTTTATTTCAC",
                "GGGATCTTACGACATAAATCG",
            ],
            "start": [0, 2, 4, 5, 8, 10, 20],
            "end": [3, 5, 7, 9, 11, 14, 25],
        }
    )


@pytest.fixture
def sequential_solution():
    return ["seq1", "seq3", "seq5", "seq7"]


@pytest.fixture
def optimal_solution():
    return ["seq1", "seq4", "seq6", "seq7"]


@pytest.mark.parametrize(
    "x,y,overlapping",
    [
        ((0, 10), (0, 10), True),
        ((0, 10), (10, 20), True),
        ((0, 10), (11, 20), False),
        ((0, 10), (4, 7), True),
    ],
)
def test_is_overlapping(x, y, overlapping):
    assert is_overlapping(x, y) == overlapping
    assert is_overlapping(y, x) == overlapping


def test_greedy_model(df, sequential_solution):
    assert greedy_model(df) == sequential_solution


@pytest.mark.skipif(not GLPK_AVAILABLE, reason="GLPK solver (glpsol) not installed")
def test_optimal_model(df, optimal_solution):
    assert optimal_model(df, 1) == optimal_solution
    assert optimal_model(df, 10) == optimal_solution


@pytest.mark.skipif(not GLPK_AVAILABLE, reason="GLPK solver (glpsol) not installed")
def test_run_optimal_as_block(df, optimal_solution):
    task = OptimizeProbeCoverage()
    task.df = df
    assigned = task.run_optimal(1, 1)
    assert assigned == optimal_solution
    assert task.df["block"].nunique() == 2


@pytest.mark.parametrize(
    "seq1,seq2,percentage,binding",
    [
        ("ATGC", "ATGC", 0.75, False),  # same seq
        ("ATGC", "ATGCA", 0.75, False),  # same but longer
        ("ATGC", "TACG", 0.75, True),  # complement
        ("ATGC", "GCAT", 0.75, True),  # reverse complement
        ("ATGCA", "TAGGT", 0.75, True),  # partially complement
        ("ATGCA", "TAGGT", 0.95, False),  # partially complement
    ],
)
def test_is_binding(seq1, seq2, percentage, binding):
    assert is_binding(seq1, seq2, percentage) == binding


@pytest.mark.parametrize(
    "invalid_percentage",
    [
        -0.5,  # negative
        -0.01,  # slightly negative
        1.01,  # slightly over 1
        1.5,  # over 1
        2.0,  # way over 1
    ],
)
def test_is_binding_invalid_percentage(invalid_percentage):
    """Test that invalid match_percentage values raise ValueError."""
    with pytest.raises(ValueError, match="between 0 and 100"):
        is_binding("ATGC", "GCAT", invalid_percentage)


def test_greedy_model_no_overlapping_assignments():
    """Test that greedy model doesn't assign overlapping probes.

    This test exposes the bug where greedy only checks the LAST assigned probe,
    not ALL assigned probes. In this case:
    - seq1 (0-10) is assigned first
    - seq2 (20-30) doesn't overlap with seq1, so it's assigned
    - seq3 (5-15) doesn't overlap with seq2 (last assigned), but DOES overlap with seq1
    - The buggy algorithm would incorrectly add seq3

    The correct result should NOT include seq3 since it overlaps with seq1.
    """
    df = pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3"],
            "sequence": ["A" * 11, "T" * 11, "G" * 11],
            "start": [0, 20, 5],
            "end": [10, 30, 15],
        }
    )
    assigned = greedy_model(df)

    # seq3 overlaps with seq1 (positions 5-10), so should NOT be in result
    assert "seq1" in assigned
    assert "seq2" in assigned
    assert "seq3" not in assigned, "seq3 overlaps with seq1 and should not be assigned"


def test_filter_binding_probes(df):
    task = OptimizeProbeCoverage()
    task.df = df
    assigned = df["name"].values

    filtered = task.filter_binding_probes(assigned, 0.75)

    assert len(filtered) <= len(assigned)
    assert all([f in assigned for f in filtered])


def test_fill_coverage_gaps():
    """Gap filler should place a probe in an uncovered region."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": ["A" * 10, "T" * 10, "G" * 10],
            "start": [0, 5, 50],
            "end": [10, 15, 60],
            "length": [10, 10, 10],
        }
    )
    # p1 and p3 assigned, gap at 10-50. p2 (5-15) doesn't fully fit in gap.
    assigned = ["p1", "p3"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    # p2 overlaps with p1, so it shouldn't be added
    assert result == ["p1", "p3"]


def test_fill_coverage_gaps_fills():
    """Gap filler should add a probe when it fits cleanly."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": ["A" * 10, "T" * 10, "G" * 10],
            "start": [0, 20, 50],
            "end": [10, 30, 60],
            "length": [10, 10, 10],
        }
    )
    # p1 and p3 assigned, gap at 10-50. p2 (20-30) fits in gap.
    assigned = ["p1", "p3"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert "p2" in result
    assert len(result) == 3


def test_fill_coverage_gaps_no_gaps():
    """No gaps means no changes."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2"],
            "sequence": ["A" * 10, "T" * 10],
            "start": [0, 11],
            "end": [10, 21],
            "length": [10, 10],
        }
    )
    assigned = ["p1", "p2"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert result == ["p1", "p2"]
