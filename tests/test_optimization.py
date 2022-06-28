import pandas as pd
import pytest

from eFISHent.optimization import OptimizeProbeCoverage
from eFISHent.optimization import is_overlapping
from eFISHent.optimization import greedy_model
from eFISHent.optimization import optimal_model


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"],
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
        ((0, 10), (4, 7), False),
    ],
)
def test_is_overlapping(x, y, overlapping):
    assert is_overlapping(x, y) == overlapping
    assert is_overlapping(y, x) == overlapping


def test_greedy_model(df, sequential_solution):
    assert greedy_model(df) == sequential_solution


def test_optimal_model(df, optimal_solution):
    assert optimal_model(df, 1) == optimal_solution
    assert optimal_model(df, 10) == optimal_solution


def test_run_optimal_as_block(df, optimal_solution):
    task = OptimizeProbeCoverage()
    task.df = df
    assigned = task.run_optimal(1, 1)
    assert assigned == optimal_solution
    assert task.df["block"].nunique() == 2
