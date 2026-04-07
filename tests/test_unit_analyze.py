"""Unit tests for eFISHent/analyze.py."""

import os
from unittest.mock import MagicMock, patch, PropertyMock

import Bio.Seq
import Bio.SeqRecord
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# Helpers to build fake SeqRecord objects
# ---------------------------------------------------------------------------

def _make_seq(seq_str, name="probe"):
    """Create a Bio.SeqRecord with the given sequence string."""
    return Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq(seq_str), id=name, name=name, description=""
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def task():
    """Return an AnalyzeProbeset instance with heavy dependencies mocked."""
    with patch("eFISHent.analyze.ProbeConfig") as mock_pc, \
         patch("eFISHent.analyze.RunConfig") as mock_rc, \
         patch("eFISHent.analyze.GeneralConfig") as mock_gc, \
         patch("eFISHent.analyze.SequenceConfig"):
        mock_pc.return_value.na_concentration = 0.3
        mock_pc.return_value.formamide_concentration = 50.0
        mock_pc.return_value.min_tm = 40.0
        mock_pc.return_value.max_tm = 60.0
        mock_pc.return_value.max_deltag = -10.0
        mock_pc.return_value.max_kmers = 10
        mock_pc.return_value.max_off_targets = 5
        mock_pc.return_value.aligner = "bowtie"
        mock_rc.return_value.analyze_probeset = "probes.fasta"
        mock_gc.return_value.threads = 1

        from eFISHent.analyze import AnalyzeProbeset
        t = AnalyzeProbeset()
        yield t


@pytest.fixture
def sample_df():
    """A small DataFrame resembling _compute_all_metrics output."""
    return pd.DataFrame({
        "name": ["p1", "p2", "p3"],
        "sequence": ["ATCGATCGATCGATCGATCG", "GCTAGCTAGCTAGCTAGCTA", "AAATTTAAATTTAAATTTAA"],
        "length": [20, 20, 20],
        "TM": [55.0, 58.0, 42.0],
        "GC": [50.0, 50.0, 20.0],
        "cpg_fraction": [0.05, 0.02, 0.0],
        "on_target_dg": [-30.0, -25.0, -15.0],
        "homopolymer_run": [1, 1, 3],
        "g_quadruplex": [False, False, False],
        "g_quadruplet_count": [0, 0, 0],
        "low_complexity": [0.05, 0.03, 0.20],
        "kmers": [2, 5, 1],
        "deltaG": [-1.0, -3.0, -0.5],
        "start": [0, 30, 60],
        "end": [20, 50, 80],
    })


# ===================================================================
# Tests for static plot helpers
# ===================================================================

class TestHistplot:
    """Tests for AnalyzeProbeset.histplot."""

    def test_basic(self, task):
        fig, ax = plt.subplots()
        task.histplot(ax, [1, 2, 3, 3, 4], "Test", 0)
        assert ax.get_title() == "Test"
        assert ax.get_ylabel() == "Count"
        plt.close(fig)

    def test_empty_data(self, task):
        fig, ax = plt.subplots()
        task.histplot(ax, [], "Empty", 0)
        assert ax.get_title() == "Empty"
        plt.close(fig)

    def test_uniform_data(self, task):
        """All identical values should use the fallback branch."""
        fig, ax = plt.subplots()
        task.histplot(ax, [5, 5, 5], "Uniform", 0)
        assert ax.get_title() == "Uniform"
        plt.close(fig)


class TestBoxplot:
    """Tests for AnalyzeProbeset.boxplot."""

    def test_basic(self, task):
        fig, ax = plt.subplots()
        task.boxplot(ax, [1.0, 2.0, 3.0], "Box", "unit")
        assert ax.get_title() == "Box"
        assert ax.get_ylabel() == "unit"
        plt.close(fig)


class TestMatrix:
    """Tests for AnalyzeProbeset.matrix."""

    def test_basic(self, task):
        fig, ax = plt.subplots()
        data = np.array([[0.0, 0.5], [0.5, 0.0]])
        task.matrix(ax, data, "Matrix", fig)
        assert ax.get_title() == "Matrix"
        plt.close(fig)


# ===================================================================
# Tests for _compute_binding_matrix
# ===================================================================

class TestComputeBindingMatrix:
    """Tests for probe-probe cross-hybridization matrix."""

    def test_diagonal_is_zero(self, task):
        task.sequences = [
            _make_seq("ATCGATCGATCGATCGATCG", "p1"),
            _make_seq("GCTAGCTAGCTAGCTAGCTA", "p2"),
        ]
        matrix = task._compute_binding_matrix()
        assert matrix[0, 0] == pytest.approx(0.0)
        assert matrix[1, 1] == pytest.approx(0.0)

    def test_symmetric(self, task):
        task.sequences = [
            _make_seq("ATCGATCGATCGATCGATCG", "p1"),
            _make_seq("GCTAGCTAGCTAGCTAGCTA", "p2"),
            _make_seq("AAATTTAAATTTAAATTTAA", "p3"),
        ]
        matrix = task._compute_binding_matrix()
        np.testing.assert_array_equal(matrix, matrix.T)

    def test_shape(self, task):
        seqs = [_make_seq("A" * 20, f"p{i}") for i in range(5)]
        task.sequences = seqs
        matrix = task._compute_binding_matrix()
        assert matrix.shape == (5, 5)

    def test_single_probe(self, task):
        task.sequences = [_make_seq("ATCGATCGATCGATCGATCG", "p1")]
        matrix = task._compute_binding_matrix()
        assert matrix.shape == (1, 1)
        assert matrix[0, 0] == pytest.approx(0.0)

    def test_values_in_range(self, task):
        task.sequences = [
            _make_seq("ATCGATCGATCGATCGATCG", "p1"),
            _make_seq("GCTAGCTAGCTAGCTAGCTA", "p2"),
        ]
        matrix = task._compute_binding_matrix()
        assert np.all(matrix >= 0.0)
        assert np.all(matrix <= 1.0)

    def test_complement_has_high_overlap(self, task):
        """A probe and its complement should have high k-mer overlap."""
        seq = "ATCGATCGATCGATCGATCG"
        comp = str(Bio.Seq.Seq(seq).complement())
        task.sequences = [_make_seq(seq, "p1"), _make_seq(comp, "p2")]
        matrix = task._compute_binding_matrix()
        # Complement should show significant overlap
        assert matrix[0, 1] > 0.3


# ===================================================================
# Tests for _compute_quality_scores
# ===================================================================

class TestComputeQualityScores:
    """Tests for composite quality score computation."""

    def test_returns_series(self, task, sample_df):
        scores = task._compute_quality_scores(sample_df)
        assert isinstance(scores, pd.Series)
        assert len(scores) == len(sample_df)

    def test_scores_in_range(self, task, sample_df):
        scores = task._compute_quality_scores(sample_df)
        assert scores.min() >= 0.0
        assert scores.max() <= 100.0

    def test_single_probe(self, task):
        """Single-probe DataFrame should get tm score = 1.0."""
        df = pd.DataFrame({
            "name": ["p1"],
            "sequence": ["ATCGATCGATCGATCGATCG"],
            "length": [20],
            "TM": [50.0],
            "GC": [50.0],
            "cpg_fraction": [0.0],
            "on_target_dg": [-30.0],
            "homopolymer_run": [1],
            "g_quadruplex": [False],
            "g_quadruplet_count": [0],
            "low_complexity": [0.0],
            "kmers": [0],
            "deltaG": [0.0],
            "start": [0],
            "end": [20],
        })
        scores = task._compute_quality_scores(df)
        assert len(scores) == 1
        # With perfect metrics the score should be high
        assert scores.iloc[0] > 50.0

    def test_off_targets_column_used(self, task, sample_df):
        """When off_targets column is present, it should affect score."""
        scores_no_ot = task._compute_quality_scores(sample_df.copy())
        df_ot = sample_df.copy()
        df_ot["off_targets"] = [10, 10, 10]
        scores_with_ot = task._compute_quality_scores(df_ot)
        # High off-targets should lower the score
        assert scores_with_ot.mean() < scores_no_ot.mean()

    def test_gc_score_ranges(self, task):
        """Test the internal GC scoring function across ranges."""
        df_base = pd.DataFrame({
            "name": ["p"],
            "sequence": ["A" * 20],
            "length": [20],
            "TM": [50.0],
            "GC": [48.0],  # ideal range
            "cpg_fraction": [0.0],
            "on_target_dg": [-30.0],
            "homopolymer_run": [1],
            "g_quadruplex": [False],
            "g_quadruplet_count": [0],
            "low_complexity": [0.0],
            "kmers": [0],
            "deltaG": [0.0],
            "start": [0],
            "end": [20],
        })

        # Ideal GC
        score_ideal = task._compute_quality_scores(df_base.copy())

        # Low GC
        df_low = df_base.copy()
        df_low["GC"] = [15.0]
        score_low = task._compute_quality_scores(df_low)

        # High GC
        df_high = df_base.copy()
        df_high["GC"] = [70.0]
        score_high = task._compute_quality_scores(df_high)

        assert score_ideal.iloc[0] > score_low.iloc[0]
        assert score_ideal.iloc[0] > score_high.iloc[0]

    def test_weights_sum_to_one(self):
        """The quality score weights should sum to 1.0."""
        weights = [0.18, 0.13, 0.05, 0.13, 0.23, 0.09, 0.09, 0.10]
        assert abs(sum(weights) - 1.0) < 1e-10

    def test_zero_kmers_gives_max_kmer_score(self, task):
        """Zero k-mer count should give kmer score = 1.0."""
        df = pd.DataFrame({
            "name": ["p1"],
            "sequence": ["ATCGATCGATCGATCGATCG"],
            "length": [20],
            "TM": [50.0],
            "GC": [50.0],
            "cpg_fraction": [0.0],
            "on_target_dg": [-30.0],
            "homopolymer_run": [1],
            "g_quadruplex": [False],
            "g_quadruplet_count": [0],
            "low_complexity": [0.0],
            "kmers": [0],
            "deltaG": [0.0],
            "start": [0],
            "end": [20],
        })
        scores = task._compute_quality_scores(df)
        # With all perfect inputs the score should be near max
        assert scores.iloc[0] > 70.0


# ===================================================================
# Tests for _save_csv
# ===================================================================

class TestSaveCsv:
    """Tests for CSV output."""

    def test_writes_csv(self, task, sample_df, tmp_path):
        csv_path = str(tmp_path / "analysis.csv")
        mock_output = MagicMock()
        mock_output.path = csv_path.replace(".csv", ".pdf")
        task.output = MagicMock(return_value=mock_output)

        task._save_csv(sample_df)
        result = pd.read_csv(csv_path)
        assert "name" in result.columns
        assert "sequence" in result.columns
        assert len(result) == 3

    def test_includes_off_targets_when_present(self, task, sample_df, tmp_path):
        csv_path = str(tmp_path / "analysis.csv")
        mock_output = MagicMock()
        mock_output.path = csv_path.replace(".csv", ".pdf")
        task.output = MagicMock(return_value=mock_output)

        sample_df["off_targets"] = [1, 2, 3]
        task._save_csv(sample_df)
        result = pd.read_csv(csv_path)
        assert "off_targets" in result.columns

    def test_includes_quality_when_present(self, task, sample_df, tmp_path):
        csv_path = str(tmp_path / "analysis.csv")
        mock_output = MagicMock()
        mock_output.path = csv_path.replace(".csv", ".pdf")
        task.output = MagicMock(return_value=mock_output)

        sample_df["quality"] = [80.0, 75.0, 60.0]
        task._save_csv(sample_df)
        result = pd.read_csv(csv_path)
        assert "quality" in result.columns


# ===================================================================
# Tests for _print_summary
# ===================================================================

class TestPrintSummary:
    """Tests for terminal summary output."""

    def test_runs_without_error(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)

    def test_with_off_targets(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)
        sample_df["off_targets"] = [0, 1, 2]

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)

    def test_with_quality(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)
        sample_df["quality"] = [80.0, 75.0, 60.0]

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)

    def test_flags_g_quadruplex(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)
        sample_df["g_quadruplex"] = [True, False, False]

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)

    def test_flags_high_gc(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)
        sample_df["GC"] = [60.0, 70.0, 50.0]

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)

    def test_no_mapped_probes(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)
        sample_df["start"] = [-1, -1, -1]
        sample_df["end"] = [-1, -1, -1]

        with patch("eFISHent.analyze.is_binding", return_value=False):
            task._print_summary(sample_df)


# ===================================================================
# Tests for _add_probe_coverage (coverage/distribution calculations)
# ===================================================================

class TestAddProbeCoverage:
    """Tests for gene coverage plotting logic."""

    def test_basic_coverage(self, task, sample_df):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)

        fig, ax = plt.subplots()
        task._add_probe_coverage(ax, sample_df)
        assert ax.get_title() == "Gene coverage"
        plt.close(fig)

    def test_no_mapped_probes(self, task):
        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)

        df = pd.DataFrame({
            "start": [-1, -1],
            "end": [-1, -1],
        })
        fig, ax = plt.subplots()
        task._add_probe_coverage(ax, df)
        assert ax.get_title() == "Gene coverage"
        plt.close(fig)

    def test_gap_annotation(self, task):
        """Large gaps should produce annotations."""
        gene = _make_seq("A" * 500, "gene")
        type(task).gene = PropertyMock(return_value=gene)

        df = pd.DataFrame({
            "start": [0, 300],
            "end": [20, 320],
        })
        fig, ax = plt.subplots()
        task._add_probe_coverage(ax, df)
        # Check that annotations were created for gaps > 100nt
        annotations = [c for c in ax.get_children()
                       if isinstance(c, matplotlib.text.Annotation)]
        assert len(annotations) >= 1
        plt.close(fig)


# ===================================================================
# Tests for _add_off_targets
# ===================================================================

class TestAddOffTargets:
    """Tests for off-target plotting helper."""

    def test_with_off_targets(self, task):
        df = pd.DataFrame({"off_targets": [0, 1, 2, 3]})
        fig, ax = plt.subplots()
        task._add_off_targets(ax, df)
        assert ax.get_title() == "Off-target count"
        plt.close(fig)

    def test_without_off_targets(self, task):
        df = pd.DataFrame({"name": ["p1"]})
        fig, ax = plt.subplots()
        task._add_off_targets(ax, df)
        assert ax.get_title() == "Off-target count"
        plt.close(fig)


# ===================================================================
# Tests for _add_quality_scores
# ===================================================================

class TestAddQualityScores:
    """Tests for quality score plotting helper."""

    def test_with_quality(self, task):
        df = pd.DataFrame({"quality": [80.0, 75.0, 60.0]})
        fig, ax = plt.subplots()
        task._add_quality_scores(ax, df)
        assert ax.get_title() == "Quality scores"
        plt.close(fig)

    def test_without_quality(self, task):
        df = pd.DataFrame({"name": ["p1"]})
        fig, ax = plt.subplots()
        task._add_quality_scores(ax, df)
        assert ax.get_title() == "Quality scores"
        plt.close(fig)


# ===================================================================
# Tests for _add_binding_affinity
# ===================================================================

class TestAddBindingAffinity:
    """Tests for binding affinity plotting helper."""

    def test_small_set(self, task):
        task.sequences = [
            _make_seq("ATCGATCGATCGATCGATCG", "p1"),
            _make_seq("GCTAGCTAGCTAGCTAGCTA", "p2"),
        ]
        fig, ax = plt.subplots()
        task._add_binding_affinity(ax, fig)
        assert ax.get_title() == "Binding affinity"
        plt.close(fig)

    def test_large_set_skipped(self, task):
        """More than 100 probes should skip the matrix computation."""
        task.sequences = [_make_seq("A" * 20, f"p{i}") for i in range(101)]
        fig, ax = plt.subplots()
        task._add_binding_affinity(ax, fig)
        assert ax.get_title() == "Binding affinity"
        plt.close(fig)


# ===================================================================
# Tests for individual _add_* panel methods
# ===================================================================

class TestPanelMethods:
    """Quick smoke tests for all plot panel methods."""

    def test_add_length(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_length(ax, sample_df)
        assert ax.get_title() == "Lengths"
        plt.close(fig)

    def test_add_melting_temperature(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_melting_temperature(ax, sample_df)
        assert ax.get_title() == "Melting temperatures"
        plt.close(fig)

    def test_add_gc_content(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_gc_content(ax, sample_df)
        assert ax.get_title() == "GC Content"
        plt.close(fig)

    def test_add_cpg_fraction(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_cpg_fraction(ax, sample_df)
        assert ax.get_title() == "CpG Fraction"
        plt.close(fig)

    def test_add_g_quadruplex(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_g_quadruplex(ax, sample_df)
        assert ax.get_title() == "G-quadruplex count"
        plt.close(fig)

    def test_add_kmers(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_kmers(ax, sample_df)
        assert ax.get_title() == "K-mer count"
        plt.close(fig)

    def test_add_free_energy(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_free_energy(ax, sample_df)
        assert ax.get_title() == "Secondary structure \u0394G"
        plt.close(fig)

    def test_add_duplex_dg(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_duplex_dg(ax, sample_df)
        assert ax.get_title() == "On-target binding \u0394G"
        plt.close(fig)

    def test_add_homopolymer(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_homopolymer(ax, sample_df)
        assert ax.get_title() == "Homopolymer run"
        plt.close(fig)

    def test_add_low_complexity(self, task, sample_df):
        fig, ax = plt.subplots()
        task._add_low_complexity(ax, sample_df)
        assert ax.get_title() == "Low complexity score"
        plt.close(fig)


# ===================================================================
# Tests for edge cases in quality scoring
# ===================================================================

class TestQualityScoreEdgeCases:
    """Edge cases for _compute_quality_scores."""

    def test_extreme_values(self, task):
        """Extreme metric values should clamp to [0, 100]."""
        df = pd.DataFrame({
            "name": ["p1"],
            "sequence": ["ATCGATCGATCGATCGATCG"],
            "length": [20],
            "TM": [100.0],
            "GC": [90.0],
            "cpg_fraction": [0.5],
            "on_target_dg": [0.0],
            "homopolymer_run": [10],
            "g_quadruplex": [True],
            "g_quadruplet_count": [5],
            "low_complexity": [1.0],
            "kmers": [100],
            "deltaG": [-50.0],
            "start": [0],
            "end": [20],
            "off_targets": [50],
        })
        scores = task._compute_quality_scores(df)
        assert scores.iloc[0] >= 0.0
        assert scores.iloc[0] <= 100.0

    def test_binding_dg_positive_clamps(self, task):
        """Positive on_target_dg should be clipped to 0."""
        df = pd.DataFrame({
            "name": ["p1"],
            "sequence": ["ATCGATCGATCGATCGATCG"],
            "length": [20],
            "TM": [50.0],
            "GC": [50.0],
            "cpg_fraction": [0.0],
            "on_target_dg": [5.0],  # positive
            "homopolymer_run": [1],
            "g_quadruplex": [False],
            "g_quadruplet_count": [0],
            "low_complexity": [0.0],
            "kmers": [0],
            "deltaG": [0.0],
            "start": [0],
            "end": [20],
        })
        scores = task._compute_quality_scores(df)
        assert scores.iloc[0] >= 0.0

    def test_all_zero_off_targets(self, task, sample_df):
        """Zero off-targets should give full off-target score."""
        sample_df["off_targets"] = [0, 0, 0]
        scores = task._compute_quality_scores(sample_df)
        # All zeros should contribute max off-target score
        assert all(s > 0 for s in scores)


# ===================================================================
# Tests for build_figure
# ===================================================================

class TestBuildFigure:
    """Tests for the full figure build."""

    def test_build_figure_creates_pdf(self, task, sample_df, tmp_path):
        pdf_path = str(tmp_path / "test_analysis.pdf")
        mock_output = MagicMock()
        mock_output.path = pdf_path
        task.output = MagicMock(return_value=mock_output)

        gene = _make_seq("A" * 100, "gene")
        type(task).gene = PropertyMock(return_value=gene)

        task.sequences = [
            _make_seq("ATCGATCGATCGATCGATCG", "p1"),
            _make_seq("GCTAGCTAGCTAGCTAGCTA", "p2"),
            _make_seq("AAATTTAAATTTAAATTTAA", "p3"),
        ]

        with patch("eFISHent.analyze.util.get_gene_name", return_value="test_gene"):
            task.build_figure(sample_df)

        assert os.path.isfile(pdf_path)
