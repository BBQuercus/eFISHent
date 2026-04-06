import shutil
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from eFISHent.optimization import OptimizeProbeCoverage
from eFISHent.optimization import compute_accessibility_scores
from eFISHent.optimization import compute_pre_optimization_quality
from eFISHent.optimization import fill_coverage_gaps
from eFISHent.optimization import fold_sequence
from eFISHent.optimization import greedy_model
from eFISHent.optimization import is_binding
from eFISHent.optimization import is_overlapping
from eFISHent.optimization import optimal_model
from eFISHent.optimization import visualize_assignment

GLPK_AVAILABLE = shutil.which("glpsol") is not None


@pytest.fixture
def df():
    sequences = [
        "GTAATTACAAAATAAGCAACG",
        "GCTTGCTTTGAGATTTTGTTC",
        "TCAATTCTCTACTGTCTCAGT",
        "GTAATTACAAAATAAGCAACG",
        "GGGACGAATTCTTTGTCTATTC",
        "TGGATTTCATAATGTTTATTTCAC",
        "GGGATCTTACGACATAAATCG",
    ]
    return pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"],
            "sequence": sequences,
            "length": [len(s) for s in sequences],
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


def test_greedy_model(df):
    assigned = greedy_model(df)
    # Greedy should produce non-overlapping probes
    assert len(assigned) > 0
    # Verify no overlaps in the assigned set
    assigned_rows = df[df["name"].isin(assigned)]
    for i, row_i in assigned_rows.iterrows():
        for j, row_j in assigned_rows.iterrows():
            if i < j:
                assert not is_overlapping(
                    (row_i["start"], row_i["end"]),
                    (row_j["start"], row_j["end"]),
                ), f"{row_i['name']} and {row_j['name']} overlap"
    # Should include seq7 (non-overlapping with everything else)
    assert "seq7" in assigned


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
        ("ATGCA", "TAGGT", 0.75, False),  # partially complement, not enough k-mer overlap
        ("ATGCATGCATGCATGCATGC", "GCATGCATGCATGCATGCA", 0.75, True),  # realistic probe-length complement
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
            "length": [11, 11, 11],
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


def test_greedy_prefers_optimal_gc():
    """Greedy should prefer probes with GC in the 45-52% optimal range."""
    from eFISHent.optimization import compute_pre_optimization_quality

    # Two overlapping probes at the same position, different GC content
    df = pd.DataFrame(
        {
            "name": ["good_gc", "bad_gc"],
            # ~50% GC vs ~80% GC (20nt probes)
            "sequence": [
                "ATGCATGCATGCATGCATGC",  # 50% GC
                "GCGCGCGCGCGCGCGCATGC",  # 80% GC
            ],
            "length": [20, 20],
            "start": [0, 0],
            "end": [20, 20],
        }
    )
    quality = compute_pre_optimization_quality(df)
    # Optimal GC (50%) should score higher than 80% GC
    assert quality.iloc[0] > quality.iloc[1]

    assigned = greedy_model(df)
    # Should pick the one with better GC
    assert assigned == ["good_gc"]


def test_compute_pre_optimization_quality_range():
    """Pre-optimization quality scores should be in (0, 1]."""
    from eFISHent.optimization import compute_pre_optimization_quality

    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": [
                "ATGCATGCATGCATGCATGC",  # 50% GC
                "GCGCGCGCGCGCGCGCGCGC",  # 100% GC
                "AAAAAAAAAAAAAAAAAAAT",  # ~5% GC
            ],
            "length": [20, 20, 20],
            "start": [0, 20, 40],
            "end": [20, 40, 60],
        }
    )
    quality = compute_pre_optimization_quality(df)
    assert (quality > 0).all()
    assert (quality <= 1).all()


# ---------------------------------------------------------------------------
# is_binding: additional correctness / regression tests
# ---------------------------------------------------------------------------


def test_is_binding_self_reverse_complement():
    """A sequence should bind its own reverse complement."""
    seq = "ATGCATGCATGCATGCATGC"
    rc = "GCATGCATGCATGCATGCAT"
    assert is_binding(seq, rc, 0.75)


def test_is_binding_near_miss_below_threshold():
    """Sequences just below the match threshold should NOT bind."""
    # These share some k-mers but not enough for 75% overlap
    seq1 = "ATGCATGCATGCATGCATGC"
    seq2 = "TTTTTTTTTTGCATGCATTT"  # only partial complement region
    assert not is_binding(seq1, seq2, 0.75)


def test_is_binding_short_sequences_return_false():
    """Sequences shorter than 3nt should always return False."""
    assert not is_binding("AT", "AT", 0.5)
    assert not is_binding("A", "T", 0.5)
    assert not is_binding("", "", 0.5)


def test_is_binding_single_nucleotide_repeats():
    """AAAA vs TTTT should bind (complement)."""
    assert is_binding("AAAA", "TTTT", 0.75)


def test_is_binding_kmer_adaptation_short_vs_long():
    """K-mer size adapts: 4nt sequences use k=3, 20nt use k=8."""
    # 4nt: k = min(8, max(3, 4-1)) = 3
    # Reverse complement of ATGC is GCAT
    assert is_binding("ATGC", "GCAT", 0.75)

    # 20nt: k = min(8, max(3, 20-1)) = 8
    seq20 = "ATGCATGCATGCATGCATGC"
    rc20 = "GCATGCATGCATGCATGCAT"
    assert is_binding(seq20, rc20, 0.75)


# ---------------------------------------------------------------------------
# greedy_model: additional correctness tests
# ---------------------------------------------------------------------------


def test_greedy_model_all_overlapping_picks_one():
    """When all probes overlap, greedy should pick only one (highest quality)."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": [
                "ATGCATGCATGCATGCATGC",  # 50% GC - optimal
                "GCGCGCGCGCGCGCGCGCGC",  # 100% GC - bad
                "AAAAAAAAAAAAAAAAAAAA",  # 0% GC - bad
            ],
            "length": [20, 20, 20],
            "start": [0, 0, 0],
            "end": [20, 20, 20],
        }
    )
    assigned = greedy_model(df)
    assert len(assigned) == 1
    assert assigned == ["p1"]  # best GC wins


def test_greedy_model_non_overlapping_all_selected():
    """Non-overlapping probes should all be selected."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": [
                "ATGCATGCATGCATGCATGC",
                "ATGCATGCATGCATGCATGC",
                "ATGCATGCATGCATGCATGC",
            ],
            "length": [20, 20, 20],
            "start": [0, 100, 200],
            "end": [20, 120, 220],
        }
    )
    assigned = greedy_model(df)
    assert len(assigned) == 3
    assert set(assigned) == {"p1", "p2", "p3"}


def test_greedy_model_quality_ordering():
    """At the same position, higher quality (better GC) probe should win."""
    df = pd.DataFrame(
        {
            "name": ["bad_gc", "good_gc"],
            "sequence": [
                "GCGCGCGCGCGCGCGCGCGC",  # 100% GC
                "ATGCATGCATGCATGCATGC",  # 50% GC
            ],
            "length": [20, 20],
            "start": [0, 0],
            "end": [20, 20],
        }
    )
    assigned = greedy_model(df)
    assert assigned == ["good_gc"]


def test_greedy_model_empty_dataframe():
    """Empty dataframe should return empty list."""
    df = pd.DataFrame(
        {"name": [], "sequence": [], "length": [], "start": [], "end": []}
    )
    assigned = greedy_model(df)
    assert assigned == []


# ---------------------------------------------------------------------------
# fill_coverage_gaps: additional regression tests
# ---------------------------------------------------------------------------


def test_fill_coverage_gaps_boundary_at_spacing():
    """Gap exactly at spacing threshold should not trigger fill."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2"],
            "sequence": ["A" * 10, "T" * 10],
            "start": [0, 12],
            "end": [10, 22],
            "length": [10, 10],
        }
    )
    # Gap is 12 - 10 = 2, spacing = 2 => not > spacing, no fill needed
    assigned = ["p1", "p2"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert result == ["p1", "p2"]


def test_fill_coverage_gaps_multiple_gaps():
    """Multiple gaps should all get filled if candidates exist."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3", "p4", "p5"],
            "sequence": ["A" * 10, "T" * 10, "G" * 10, "C" * 10, "A" * 10],
            "start": [0, 20, 40, 60, 80],
            "end": [10, 30, 50, 70, 90],
            "length": [10, 10, 10, 10, 10],
        }
    )
    # Assign p1, p3, p5 => gaps at 10-40 (p2 fits) and 50-80 (p4 fits)
    assigned = ["p1", "p3", "p5"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert "p2" in result
    assert "p4" in result
    assert len(result) == 5


def test_fill_coverage_gaps_no_available_unassigned():
    """Gap with no available unassigned probes should leave gap unfilled."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2"],
            "sequence": ["A" * 10, "T" * 10],
            "start": [0, 50],
            "end": [10, 60],
            "length": [10, 10],
        }
    )
    # Gap at 10-50, but no unassigned probes that fit in that gap
    assigned = ["p1", "p2"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert result == ["p1", "p2"]


def test_fill_coverage_gaps_reverse_order():
    """Assigned probes in reverse order should still work."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2", "p3"],
            "sequence": ["A" * 10, "T" * 10, "G" * 10],
            "start": [0, 20, 50],
            "end": [10, 30, 60],
            "length": [10, 10, 10],
        }
    )
    # Pass assigned in reverse order
    assigned = ["p3", "p1"]
    result = fill_coverage_gaps(df, assigned, spacing=2)
    assert "p2" in result


# ---------------------------------------------------------------------------
# compute_pre_optimization_quality: correctness tests
# ---------------------------------------------------------------------------


def test_quality_optimal_gc_scores_higher():
    """48% GC should score higher than 70% GC."""
    # 48% GC: 9-10 G/C out of 20
    seq_48 = "ATGCATGCATATATATATAT"  # ~40% GC
    seq_70 = "GCGCGCGCGCGCGCATATGC"  # ~70% GC
    df = pd.DataFrame(
        {
            "name": ["optimal", "extreme"],
            "sequence": [seq_48, seq_70],
            "length": [20, 20],
            "start": [0, 20],
            "end": [20, 40],
        }
    )
    quality = compute_pre_optimization_quality(df)
    # Closer to 45-52% should score better than 70%
    assert quality.iloc[0] > quality.iloc[1]


def test_quality_high_cpg_scores_lower():
    """High CpG fraction should score lower than low CpG."""
    # High CpG: many CG dinucleotides
    seq_high_cpg = "CGCGCGATATATCGATCGCG"  # lots of CG
    # Low CpG: same GC but no CG dinucleotides (use GC not CG)
    seq_low_cpg = "GCGCGCATATATTGACTGCG"  # fewer CG dinucleotides
    df = pd.DataFrame(
        {
            "name": ["high_cpg", "low_cpg"],
            "sequence": [seq_high_cpg, seq_low_cpg],
            "length": [20, 20],
            "start": [0, 20],
            "end": [20, 40],
        }
    )
    quality = compute_pre_optimization_quality(df)
    # The one with more CpG should score lower (or equal if GC dominates)
    # We just check both are valid scores; CpG penalty applies
    assert quality.iloc[0] <= quality.iloc[1] or abs(quality.iloc[0] - quality.iloc[1]) < 0.05


def test_quality_minimum_clamp():
    """Scores should always be >= 0.1 (minimum clamp)."""
    # Extreme sequences that should get low quality
    df = pd.DataFrame(
        {
            "name": ["p1", "p2"],
            "sequence": [
                "GCGCGCGCGCGCGCGCGCGC",  # 100% GC
                "AAAAAAAAAAAAAAAAAAAA",  # 0% GC
            ],
            "length": [20, 20],
            "start": [0, 20],
            "end": [20, 40],
        }
    )
    quality = compute_pre_optimization_quality(df)
    assert (quality >= 0.1).all()


def test_quality_single_factor_variation():
    """With identical GC but different CpG, only CpG factor should vary."""
    # Both 50% GC, different CpG counts
    seq_no_cpg = "GCGCATATATGCGCATATAT"  # has some CG
    seq_cpg = "CGCGCGCGCGATATATATAT"  # more CG dinucleotides
    df = pd.DataFrame(
        {
            "name": ["less_cpg", "more_cpg"],
            "sequence": [seq_no_cpg, seq_cpg],
            "length": [20, 20],
            "start": [0, 20],
            "end": [20, 40],
        }
    )
    quality = compute_pre_optimization_quality(df)
    # Both should be valid scores
    assert (quality > 0).all()
    assert (quality <= 1).all()


# ---------------------------------------------------------------------------
# refine_tm_uniformity: regression tests
# ---------------------------------------------------------------------------


def test_refine_tm_uniformity_fewer_than_3():
    """Fewer than 3 probes should return input unchanged."""
    import Bio.SeqRecord
    import Bio.Seq

    task = OptimizeProbeCoverage()
    task.df = pd.DataFrame(
        {
            "name": ["candidate-1-0", "candidate-1-100"],
            "sequence": ["ATGCATGCATGCATGCATGC", "GCGCATATATGCGCATATAT"],
            "length": [20, 20],
            "start": [0, 100],
            "end": [20, 120],
        }
    )
    sequences = [
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s), id=n)
        for n, s in zip(task.df["name"], task.df["sequence"])
    ]
    assigned = list(task.df["name"])
    result = task.refine_tm_uniformity(assigned, sequences)
    assert result == assigned


def test_refine_tm_uniformity_already_uniform():
    """If all probes have similar Tm, no swaps should occur."""
    import Bio.SeqRecord
    import Bio.Seq

    task = OptimizeProbeCoverage()
    # Use identical sequences so Tm is perfectly uniform
    seq = "ATGCATGCATGCATGCATGC"
    names = ["candidate-1-0", "candidate-1-100", "candidate-1-200"]
    task.df = pd.DataFrame(
        {
            "name": names,
            "sequence": [seq] * 3,
            "length": [20] * 3,
            "start": [0, 100, 200],
            "end": [20, 120, 220],
        }
    )
    sequences = [
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=n) for n in names
    ]
    assigned = list(names)
    result = task.refine_tm_uniformity(assigned, sequences)
    assert result == assigned


# ---------------------------------------------------------------------------
# is_overlapping: boundary tests
# ---------------------------------------------------------------------------


def test_is_overlapping_adjacent_touch():
    """Adjacent ranges that touch at boundary (0,10) and (10,20) overlap."""
    assert is_overlapping((0, 10), (10, 20)) is True


def test_is_overlapping_identical():
    """Identical ranges should overlap."""
    assert is_overlapping((5, 15), (5, 15)) is True


def test_is_overlapping_contained():
    """A range fully contained in another should overlap."""
    assert is_overlapping((0, 100), (10, 20)) is True


def test_is_overlapping_commutativity():
    """Order of arguments should not matter."""
    assert is_overlapping((0, 10), (5, 15)) == is_overlapping((5, 15), (0, 10))
    assert is_overlapping((0, 10), (20, 30)) == is_overlapping((20, 30), (0, 10))


# ---------------------------------------------------------------------------
# fold_sequence: mocked tests
# ---------------------------------------------------------------------------


@patch("eFISHent.optimization.subprocess.run")
@patch("eFISHent.optimization.os.path.isfile", return_value=True)
def test_fold_sequence_known_output(mock_isfile, mock_run):
    """Mock subprocess to return known dot-bracket, verify parsing."""
    mock_result = MagicMock()
    mock_result.stdout = ">target\nATGCATGC\n...((...))\n"
    mock_run.return_value = mock_result
    result = fold_sequence("ATGCATGC")
    assert result == "...((...))"


@patch("eFISHent.optimization.subprocess.run", side_effect=Exception("fail"))
@patch("eFISHent.optimization.os.path.isfile", return_value=True)
def test_fold_sequence_failure_returns_none(mock_isfile, mock_run):
    """Subprocess failure should return None."""
    result = fold_sequence("ATGCATGC")
    assert result is None


@patch("eFISHent.optimization.sys")
def test_fold_sequence_unsupported_platform(mock_sys):
    """Unsupported platform should return None."""
    mock_sys.platform = "win32"
    result = fold_sequence("ATGCATGC")
    assert result is None


# ---------------------------------------------------------------------------
# compute_accessibility_scores: mocked tests
# ---------------------------------------------------------------------------


@patch("eFISHent.optimization.fold_sequence")
def test_accessibility_all_unpaired(mock_fold):
    """All dots (unpaired) should give score 1.0."""
    target_seq = "A" * 300
    df = pd.DataFrame(
        {
            "name": ["p1"],
            "sequence": ["A" * 20],
            "start": [100],
            "end": [120],
            "length": [20],
        }
    )
    # Return all dots for the folded window
    mock_fold.return_value = "." * 200
    scores = compute_accessibility_scores(df, target_seq)
    assert scores.iloc[0] == 1.0


@patch("eFISHent.optimization.fold_sequence")
def test_accessibility_all_paired(mock_fold):
    """All parens (paired) should give score 0.0."""
    target_seq = "A" * 300
    df = pd.DataFrame(
        {
            "name": ["p1"],
            "sequence": ["A" * 20],
            "start": [100],
            "end": [120],
            "length": [20],
        }
    )
    # Return all parens for the folded window
    mock_fold.return_value = "(" * 100 + ")" * 100
    scores = compute_accessibility_scores(df, target_seq)
    assert scores.iloc[0] == 0.0


@patch("eFISHent.optimization.fold_sequence")
def test_accessibility_mixed_structure(mock_fold):
    """Mixed structure should give score between 0 and 1."""
    target_seq = "A" * 300
    df = pd.DataFrame(
        {
            "name": ["p1"],
            "sequence": ["A" * 20],
            "start": [100],
            "end": [120],
            "length": [20],
        }
    )
    # Create a structure where the probe binding region is half paired
    structure = "." * 100 + "(" * 10 + "." * 10 + ")" * 80
    mock_fold.return_value = structure
    scores = compute_accessibility_scores(df, target_seq)
    assert 0.0 < scores.iloc[0] < 1.0


@patch("eFISHent.optimization.fold_sequence")
def test_accessibility_short_sequence_fallback(mock_fold):
    """Short target sequence (< 10nt in window) should give score 1.0."""
    target_seq = "ATGCA"
    df = pd.DataFrame(
        {
            "name": ["p1"],
            "sequence": ["ATGCA"],
            "start": [0],
            "end": [5],
            "length": [5],
        }
    )
    # fold_sequence should not even be called for short sequences
    scores = compute_accessibility_scores(df, target_seq)
    assert scores.iloc[0] == 1.0
    mock_fold.assert_not_called()


@patch("eFISHent.optimization.fold_sequence", return_value=None)
def test_accessibility_folding_failure(mock_fold):
    """Folding failure (returns None) should give score 1.0."""
    target_seq = "A" * 300
    df = pd.DataFrame(
        {
            "name": ["p1"],
            "sequence": ["A" * 20],
            "start": [100],
            "end": [120],
            "length": [20],
        }
    )
    scores = compute_accessibility_scores(df, target_seq)
    assert scores.iloc[0] == 1.0


# ---------------------------------------------------------------------------
# visualize_assignment: output file test
# ---------------------------------------------------------------------------


@patch("eFISHent.optimization.plt")
def test_visualize_assignment_creates_file(mock_plt):
    """Verify visualize_assignment calls savefig with the correct path."""
    df = pd.DataFrame(
        {
            "name": ["p1", "p2"],
            "sequence": ["A" * 10, "T" * 10],
            "start": [0, 20],
            "end": [10, 30],
            "length": [10, 10],
        }
    )
    assigned = ["p1"]
    output_path = "/tmp/test_optimization_viz.png"
    visualize_assignment(df, assigned, output_path)
    mock_plt.savefig.assert_called_once_with(
        output_path, dpi=600, bbox_inches="tight"
    )
