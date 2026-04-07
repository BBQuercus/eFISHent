import os
import shutil
import tempfile
from unittest.mock import MagicMock, patch

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import luigi
import pandas as pd
import pytest

from eFISHent.cleanup import CleanUpOutput

JELLYFISH_AVAILABLE = shutil.which("jellyfish") is not None


@pytest.fixture
def task_cleanup():
    return CleanUpOutput()


# ---------------------------------------------------------------------------
# prettify_table
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not JELLYFISH_AVAILABLE, reason="jellyfish not installed")
def test_prettify_table(task_cleanup):
    class Config(luigi.Config):
        na_concentration = luigi.FloatParameter(390)
        formamide_concentration = luigi.FloatParameter(10)

    sequences = [
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATGCATGCATTAGCGCATGCATGC"), id="candidate-1-4"
        ),
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATGCATGCATGCATGCATGCATGC"), id="candidate-0-1"
        ),
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATGCATGCATTTTTTTATGCATGC"), id="candidate-6-5"
        ),
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATGCATGGGGGGGGGCATGCATGC"), id="candidate-8-20"
        ),
    ]

    # Create temp alignment CSV file with qname column
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("qname,flag,rname,pos\n")
        f.write("candidate-1-4,0,chr1,100\n")
        f.write("candidate-0-1,0,chr1,200\n")
        f.write("candidate-6-5,0,chr1,300\n")
        f.write("candidate-8-20,0,chr1,400\n")
        alignment_path = f.name

    try:
        output = task_cleanup.prettify_table(
            sequences,
            basename="test",
            jellyfish_path="./tests/data/sacCer3_15.jf",
            alignment_path=alignment_path,
            config=Config,
        )
        columns = [
            "name",
            "sequence",
            "length",
            "start",
            "end",
            "GC",
            "TM",
            "deltaG",
            "kmers",
            "count",
        ]
        for col in columns:
            assert col in output.columns
        assert len(output) == len(sequences)
        assert output.loc[0, "name"] == "test-1"
        assert output.loc[len(sequences) - 1, "name"] == f"test-{len(sequences)}"
        assert output.loc[0, "GC"] == 50
        assert output.loc[3, "GC"] == pytest.approx(66.67)
    finally:
        os.unlink(alignment_path)


def test_prettify_sequences(task_cleanup):
    data = {
        "name": {0: "test-1", 1: "test-2", 2: "test-3", 3: "test-4"},
        "sequence": {
            0: "ATGCATGCATGCATGCATGCATGC",
            1: "ATGCATGCAAGCATAGACACATGC",
            2: "TACATGCTAGACATGCATGCATGC",
            3: "ATGCAGATCGATGGGCATGCATGC",
        },
        "length": {0: 24, 1: 24, 2: 24, 3: 24},
        "start": {0: 20, 1: 20, 2: 20, 3: 20},
    }
    df = pd.DataFrame(data)
    output = task_cleanup.prettify_sequences(df)
    assert len(output) == len(df)
    for idx, row in df.iterrows():
        assert output[idx].seq == row["sequence"]
        assert output[idx].id == row["name"]


def test_cleanup_removes_blast_flagged_probes(monkeypatch, task_cleanup):
    """Final BLAST verification should remove flagged probes from outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = os.path.join(tmpdir, "optimal.fasta")
        output_fasta = os.path.join(tmpdir, "final.fasta")
        output_csv = os.path.join(tmpdir, "final.csv")
        output_cfg = os.path.join(tmpdir, "final.txt")

        Bio.SeqIO.write(
            [
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGCATGCATGCATGCATGC"), id="candidate-0-1"),
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("TTTTCCCCAAAAGGGGTTTT"), id="candidate-1-30"),
            ],
            input_fasta,
            "fasta",
        )

        df = pd.DataFrame(
            {
                "name": ["test-1", "test-2"],
                "sequence": ["ATGCATGCATGCATGCATGC", "TTTTCCCCAAAAGGGGTTTT"],
                "length": [20, 20],
                "start": [1, 30],
                "end": [21, 50],
            }
        )

        monkeypatch.setattr(
            task_cleanup,
            "input",
            lambda: {
                "optimize": {"probes": luigi.LocalTarget(input_fasta)},
                "jellyfish": luigi.LocalTarget("/tmp/fake.jf"),
                "alignment": {"table": luigi.LocalTarget("/tmp/fake.csv")},
            },
        )
        monkeypatch.setattr(
            task_cleanup,
            "output",
            lambda: {
                "fasta": luigi.LocalTarget(output_fasta),
                "table": luigi.LocalTarget(output_csv),
                "config": luigi.LocalTarget(output_cfg),
            },
        )
        monkeypatch.setattr(task_cleanup, "prettify_table", lambda *args, **kwargs: df.copy())
        monkeypatch.setattr(task_cleanup, "prettify_configuration", lambda: "config")
        monkeypatch.setattr(task_cleanup, "_compute_summary", lambda table: {"probe_count": len(table)})
        monkeypatch.setattr(task_cleanup, "remove_intermediates", lambda: None)
        monkeypatch.setattr("eFISHent.cleanup.util.get_gene_name", lambda hashed=False: "test")
        monkeypatch.setattr(
            task_cleanup,
            "_run_blast_verification",
            lambda path: {
                "total": 2,
                "clean": 1,
                "flagged": {"test-2": 2},
                "max_expected": 0,
            },
        )
        monkeypatch.setattr("eFISHent.cleanup.util.log_stage_start", lambda *args, **kwargs: None)
        monkeypatch.setattr("eFISHent.cleanup.util.log_and_check_candidates", lambda *args, **kwargs: None)

        task_cleanup.run()

        out_df = pd.read_csv(output_csv)
        out_records = list(Bio.SeqIO.parse(output_fasta, "fasta"))

        assert out_df["name"].tolist() == ["test-1"]
        assert [record.id for record in out_records] == ["test-1"]
        assert task_cleanup._verification["flagged"] == {}


# ---------------------------------------------------------------------------
# _compute_quality_scores
# ---------------------------------------------------------------------------


class _FakeProbeConfig:
    """Minimal stand-in for ProbeConfig defaults used by quality scoring."""
    min_tm = 40.0
    max_tm = 60.0
    max_deltag = -10.0
    max_kmers = 5
    max_off_targets = 0


def _make_quality_df(**overrides):
    """Build a minimal DataFrame suitable for _compute_quality_scores."""
    base = {
        "TM": [50.0],
        "GC": [48.0],
        "cpg_fraction": [0.02],
        "deltaG": [-2.0],
        "kmers": [1],
        "count": [0],
        "txome_off_targets": [0],
        "low_complexity": [0.0],
        "accessibility": [1.0],
        "on_target_dg": [-30.0],
    }
    base.update(overrides)
    return pd.DataFrame(base)


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_perfect_probe(mock_cfg):
    """A near-ideal probe should score close to 100."""
    df = _make_quality_df()
    scores = CleanUpOutput._compute_quality_scores(df)
    assert len(scores) == 1
    assert scores.iloc[0] > 80


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_high_gc_penalty(mock_cfg):
    """GC > 55% should reduce the quality score."""
    df_good = _make_quality_df(GC=[48.0])
    df_bad = _make_quality_df(GC=[65.0])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_bad = CleanUpOutput._compute_quality_scores(df_bad).iloc[0]
    assert q_good > q_bad


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_low_gc_penalty(mock_cfg):
    """GC < 45% should also reduce the quality score."""
    df_good = _make_quality_df(GC=[48.0])
    df_low = _make_quality_df(GC=[25.0])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_low = CleanUpOutput._compute_quality_scores(df_low).iloc[0]
    assert q_good > q_low


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_high_off_targets(mock_cfg):
    """Many off-targets should reduce quality."""
    df_good = _make_quality_df(count=[0])
    df_bad = _make_quality_df(count=[5])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_bad = CleanUpOutput._compute_quality_scores(df_bad).iloc[0]
    assert q_good > q_bad


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_high_kmers_penalty(mock_cfg):
    """High k-mer counts should reduce quality."""
    df_good = _make_quality_df(kmers=[0])
    df_bad = _make_quality_df(kmers=[5])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_bad = CleanUpOutput._compute_quality_scores(df_bad).iloc[0]
    assert q_good > q_bad


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_missing_optional_columns(mock_cfg):
    """Score should still work when optional columns are absent."""
    df = pd.DataFrame({
        "TM": [50.0],
        "GC": [48.0],
        "deltaG": [-2.0],
        "kmers": [1],
        "count": [0],
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    assert len(scores) == 1
    assert 0 <= scores.iloc[0] <= 100


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_multiple_probes(mock_cfg):
    """Score computation works for multiple probes at once."""
    df = pd.DataFrame({
        "TM": [50.0, 55.0, 45.0],
        "GC": [48.0, 60.0, 30.0],
        "cpg_fraction": [0.02, 0.08, 0.0],
        "deltaG": [-2.0, -8.0, 0.0],
        "kmers": [1, 4, 0],
        "count": [0, 2, 0],
        "txome_off_targets": [0, 3, 0],
        "low_complexity": [0.0, 0.2, 0.0],
        "accessibility": [1.0, 1.0, 1.0],
        "on_target_dg": [-30.0, -30.0, -30.0],
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    assert len(scores) == 3
    # First probe (ideal) should beat the second (bad gc, high off-targets)
    assert scores.iloc[0] > scores.iloc[1]


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_low_complexity_penalty(mock_cfg):
    """High low-complexity fraction should reduce quality."""
    df_good = _make_quality_df(low_complexity=[0.0])
    df_bad = _make_quality_df(low_complexity=[0.3])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_bad = CleanUpOutput._compute_quality_scores(df_bad).iloc[0]
    assert q_good > q_bad


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_cpg_penalty(mock_cfg):
    """High CpG fraction should reduce quality."""
    df_good = _make_quality_df(cpg_fraction=[0.0])
    df_bad = _make_quality_df(cpg_fraction=[0.10])
    q_good = CleanUpOutput._compute_quality_scores(df_good).iloc[0]
    q_bad = CleanUpOutput._compute_quality_scores(df_bad).iloc[0]
    assert q_good > q_bad


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_empty_dataframe(mock_cfg):
    """Empty DataFrame should return empty Series without error."""
    df = pd.DataFrame({
        "TM": pd.Series(dtype=float),
        "GC": pd.Series(dtype=float),
        "deltaG": pd.Series(dtype=float),
        "kmers": pd.Series(dtype=int),
        "count": pd.Series(dtype=int),
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    assert len(scores) == 0


# ---------------------------------------------------------------------------
# _compute_recommendation
# ---------------------------------------------------------------------------


def test_recommendation_pass():
    row = pd.Series({"quality": 85, "txome_off_targets": 0, "expression_risk": ""})
    assert CleanUpOutput._compute_recommendation(row) == "PASS"


def test_recommendation_fail_low_quality():
    row = pd.Series({"quality": 20, "txome_off_targets": 0, "expression_risk": ""})
    assert CleanUpOutput._compute_recommendation(row) == "FAIL"


def test_recommendation_flag_low_quality():
    row = pd.Series({"quality": 50, "txome_off_targets": 0, "expression_risk": ""})
    result = CleanUpOutput._compute_recommendation(row)
    assert result.startswith("FLAG")
    assert "low_quality" in result


def test_recommendation_flag_many_off_targets():
    row = pd.Series({"quality": 80, "txome_off_targets": 5, "expression_risk": ""})
    result = CleanUpOutput._compute_recommendation(row)
    assert result.startswith("FLAG")
    assert "many_off_targets" in result


def test_recommendation_flag_high_expression():
    row = pd.Series({"quality": 80, "txome_off_targets": 0, "expression_risk": "GENE:HIGH(100)"})
    result = CleanUpOutput._compute_recommendation(row)
    assert result.startswith("FLAG")
    assert "high_expression_risk" in result


def test_recommendation_flag_multiple_reasons():
    row = pd.Series({"quality": 50, "txome_off_targets": 5, "expression_risk": "GENE:HIGH(100)"})
    result = CleanUpOutput._compute_recommendation(row)
    assert result.startswith("FLAG")
    assert "low_quality" in result
    assert "many_off_targets" in result
    assert "high_expression_risk" in result


def test_recommendation_missing_keys():
    """Should gracefully handle missing keys via .get() defaults."""
    row = pd.Series({"quality": 85})
    assert CleanUpOutput._compute_recommendation(row) == "PASS"


# ---------------------------------------------------------------------------
# _compute_low_complexity_score
# ---------------------------------------------------------------------------


def test_low_complexity_homopolymer():
    """A homopolymer run has very low entropy -> high complexity score."""
    seq = Bio.Seq.Seq("AAAAAAAAAAAAAAAAAAAAAAAAA")
    score = CleanUpOutput._compute_low_complexity_score(seq)
    assert score > 0.5


def test_low_complexity_diverse():
    """A diverse sequence should have low complexity score."""
    seq = Bio.Seq.Seq("ATGCATGCATGCATGCATGCATGC")
    score = CleanUpOutput._compute_low_complexity_score(seq)
    assert score == pytest.approx(0.0)


def test_low_complexity_short_sequence():
    """Sequences shorter than the window should return 0.0."""
    seq = Bio.Seq.Seq("ATGC")
    score = CleanUpOutput._compute_low_complexity_score(seq)
    assert score == pytest.approx(0.0)


def test_low_complexity_returns_fraction():
    """Score should be between 0 and 1."""
    seq = Bio.Seq.Seq("ATGCAAAAAAAAAAATGCATGCATGC")
    score = CleanUpOutput._compute_low_complexity_score(seq)
    assert 0.0 <= score <= 1.0


# ---------------------------------------------------------------------------
# _detect_off_target_clustering
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_clustering_no_column(mock_gene, task_cleanup):
    """No off_target_genes column -> no warnings."""
    df = pd.DataFrame({"name": ["p1"], "quality": [90]})
    assert task_cleanup._detect_off_target_clustering(df) == []


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_clustering_no_hits(mock_gene, task_cleanup):
    """Empty off_target_genes -> no warnings."""
    df = pd.DataFrame({
        "name": ["p1", "p2"],
        "off_target_genes": ["", ""],
        "quality": [90, 85],
    })
    assert task_cleanup._detect_off_target_clustering(df) == []


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_clustering_below_threshold(mock_gene, task_cleanup):
    """Fewer than 5 probes hitting the same gene -> no warnings."""
    df = pd.DataFrame({
        "name": [f"p{i}" for i in range(4)],
        "off_target_genes": ["GAPDH(1)"] * 4,
        "quality": [80] * 4,
    })
    assert task_cleanup._detect_off_target_clustering(df) == []


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_clustering_detected(mock_gene, task_cleanup):
    """5+ probes hitting the same gene triggers a warning."""
    df = pd.DataFrame({
        "name": [f"p{i}" for i in range(6)],
        "off_target_genes": ["GAPDH(2)"] * 6,
        "quality": [80] * 6,
    })
    warnings = task_cleanup._detect_off_target_clustering(df)
    assert len(warnings) == 1
    assert "GAPDH" in warnings[0]


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_clustering_paralog(mock_gene, task_cleanup):
    """Paralog detection when target and off-target share a family prefix."""
    df = pd.DataFrame({
        "name": [f"p{i}" for i in range(6)],
        "off_target_genes": ["ACTG1(3)"] * 6,
        "quality": [80] * 6,
    })
    warnings = task_cleanup._detect_off_target_clustering(df)
    assert len(warnings) == 1
    assert "paralog" in warnings[0].lower()


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
def test_off_target_undruggable(mock_gene, task_cleanup):
    """>=50% of probes hitting the same off-target triggers UNDRUGGABLE."""
    n_probes = 10
    df = pd.DataFrame({
        "name": [f"p{i}" for i in range(n_probes)],
        "off_target_genes": ["ACTG1(3)"] * n_probes,
        "quality": [80] * n_probes,
    })
    warnings = task_cleanup._detect_off_target_clustering(df)
    assert len(warnings) == 1
    assert "UNDRUGGABLE" in warnings[0]
    assert "paralog" in warnings[0].lower()


@patch("eFISHent.cleanup.util.get_gene_name", return_value="MYH1")
def test_off_target_undruggable_non_paralog(mock_gene, task_cleanup):
    """>=50% of probes hitting a non-paralog off-target."""
    n_probes = 10
    df = pd.DataFrame({
        "name": [f"p{i}" for i in range(n_probes)],
        "off_target_genes": ["GAPDH(3)"] * n_probes,
        "quality": [80] * n_probes,
    })
    warnings = task_cleanup._detect_off_target_clustering(df)
    assert len(warnings) == 1
    assert "UNDRUGGABLE" in warnings[0]
    assert "paralog" not in warnings[0].lower()


# ---------------------------------------------------------------------------
# _apply_off_target_cap
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.ProbeConfig")
def test_off_target_cap_disabled(mock_cfg, task_cleanup):
    """Cap <= 0 should return the dataframe unchanged."""
    mock_cfg.return_value.max_probes_per_off_target = 0
    df = pd.DataFrame({
        "name": ["p1", "p2"],
        "off_target_genes": ["GAPDH(1)", "GAPDH(1)"],
        "quality": [80, 90],
    })
    result = task_cleanup._apply_off_target_cap(df)
    assert len(result) == 2


@patch("eFISHent.cleanup.ProbeConfig")
def test_off_target_cap_removes_lowest_quality(mock_cfg, task_cleanup):
    """When cap=1, should keep the higher quality probe."""
    mock_cfg.return_value.max_probes_per_off_target = 1
    df = pd.DataFrame({
        "name": ["p1", "p2", "p3"],
        "off_target_genes": ["GAPDH(1)", "GAPDH(1)", ""],
        "quality": [90, 50, 80],
    })
    result = task_cleanup._apply_off_target_cap(df)
    # p2 (quality=50) should be removed; p1 and p3 remain
    assert len(result) == 2
    assert "p1" in result["name"].values
    assert "p3" in result["name"].values


@patch("eFISHent.cleanup.ProbeConfig")
def test_off_target_cap_no_off_target_column(mock_cfg, task_cleanup):
    """Missing off_target_genes column returns unchanged."""
    mock_cfg.return_value.max_probes_per_off_target = 1
    df = pd.DataFrame({"name": ["p1"], "quality": [80]})
    result = task_cleanup._apply_off_target_cap(df)
    assert len(result) == 1


# ---------------------------------------------------------------------------
# prettify_section / prettify_configuration
# ---------------------------------------------------------------------------


def test_prettify_section(task_cleanup):
    """prettify_section formats a config section nicely."""
    task_cleanup.config = {
        "ProbeConfig": {"min_tm": "40.0", "max_tm": "60.0", "empty_val": "", "quoted": '""'}
    }
    result = task_cleanup.prettify_section("ProbeConfig")
    assert "Probe configuration:" in result
    assert "min-tm: 40.0" in result
    assert "max-tm: 60.0" in result
    # Empty and quoted-empty values should be excluded
    assert "empty-val" not in result
    assert "quoted" not in result


def test_prettify_configuration(task_cleanup, monkeypatch):
    """prettify_configuration merges all sections and adds the command."""
    fake_config = MagicMock()
    fake_config.sections.return_value = ["ProbeConfig", "RunConfig"]
    fake_config.__getitem__ = lambda self, key: {
        "ProbeConfig": {"min_tm": "40.0"},
        "RunConfig": {"threads": "4"},
    }[key]
    monkeypatch.setattr(luigi.configuration, "get_config", lambda: fake_config)
    monkeypatch.setattr("sys.argv", ["efishent", "--gene", "ACTB"])

    result = task_cleanup.prettify_configuration()
    assert "Probe configuration:" in result
    assert "Run configuration:" in result
    assert "Command:" in result
    assert "--gene ACTB" in result


# ---------------------------------------------------------------------------
# remove_intermediates
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.util.get_output_dir")
@patch("eFISHent.cleanup.util.get_gene_name")
def test_remove_intermediates(mock_gene, mock_outdir, task_cleanup):
    """Should remove intermediate files but keep the entrez fasta."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_outdir.return_value = tmpdir
        mock_gene.return_value = "test123"

        # Create intermediate files
        inter1 = os.path.join(tmpdir, "test123_probes.fasta")
        inter2 = os.path.join(tmpdir, "test123_alignment.csv")
        entrez = os.path.join(tmpdir, "test123_entrez.fasta")
        for f in [inter1, inter2, entrez]:
            open(f, "w").close()

        task_cleanup.remove_intermediates()

        assert not os.path.exists(inter1)
        assert not os.path.exists(inter2)
        assert os.path.exists(entrez)


@patch("eFISHent.cleanup.util.get_output_dir")
@patch("eFISHent.cleanup.util.get_gene_name")
def test_remove_intermediates_no_entrez(mock_gene, mock_outdir, task_cleanup):
    """Should not crash when entrez file does not exist."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_outdir.return_value = tmpdir
        mock_gene.return_value = "test123"

        inter1 = os.path.join(tmpdir, "test123_probes.fasta")
        open(inter1, "w").close()

        task_cleanup.remove_intermediates()
        assert not os.path.exists(inter1)


# ---------------------------------------------------------------------------
# _annotate_expression_risk
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.GeneralConfig")
@patch("eFISHent.cleanup.ProbeConfig")
def test_annotate_expression_risk_no_annotation(mock_pcfg, mock_gcfg, task_cleanup):
    """Without annotation, expression_risk column should be empty strings."""
    mock_gcfg.return_value.reference_annotation = ""
    df = pd.DataFrame({
        "name": ["p1"],
        "off_target_genes": ["GAPDH(1)"],
    })
    task_cleanup._annotate_expression_risk(df)
    assert df["expression_risk"].iloc[0] == ""


@patch("eFISHent.cleanup.GeneralConfig")
@patch("eFISHent.cleanup.ProbeConfig")
@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
@patch("eFISHent.cleanup.util.get_gene_name", return_value="TEST")
@patch("eFISHent.cleanup.build_transcript_gene_map", return_value={"ENS001": "GAPDH"})
def test_annotate_expression_risk_with_count_table(
    mock_map, mock_gene, mock_outdir, mock_pcfg, mock_gcfg, task_cleanup
):
    """With an expression table, probes with highly expressed off-targets get HIGH risk."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("gene_id,expression\n")
        f.write("ENS001,100.0\n")
        count_table = f.name

    try:
        mock_gcfg.return_value.reference_annotation = "/fake/annotation.gtf"
        mock_pcfg.return_value.encode_count_table = count_table

        df = pd.DataFrame({
            "name": ["p1"],
            "off_target_genes": ["GAPDH(2)"],
        })
        task_cleanup._annotate_expression_risk(df)
        assert "HIGH" in df["expression_risk"].iloc[0]
        assert "GAPDH" in df["expression_risk"].iloc[0]
    finally:
        os.unlink(count_table)


@patch("eFISHent.cleanup.GeneralConfig")
@patch("eFISHent.cleanup.ProbeConfig")
@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
@patch("eFISHent.cleanup.util.get_gene_name", return_value="TEST")
@patch("eFISHent.cleanup.build_transcript_gene_map", return_value={"ENS001": "GAPDH"})
def test_annotate_expression_risk_moderate(
    mock_map, mock_gene, mock_outdir, mock_pcfg, mock_gcfg, task_cleanup
):
    """Moderate expression (10-50 TPM) should get moderate risk."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("gene_id,expression\n")
        f.write("ENS001,25.0\n")
        count_table = f.name

    try:
        mock_gcfg.return_value.reference_annotation = "/fake/annotation.gtf"
        mock_pcfg.return_value.encode_count_table = count_table

        df = pd.DataFrame({
            "name": ["p1"],
            "off_target_genes": ["GAPDH(2)"],
        })
        task_cleanup._annotate_expression_risk(df)
        assert "moderate" in df["expression_risk"].iloc[0]
    finally:
        os.unlink(count_table)


@patch("eFISHent.cleanup.GeneralConfig")
@patch("eFISHent.cleanup.ProbeConfig")
@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
@patch("eFISHent.cleanup.util.get_gene_name", return_value="TEST")
@patch("eFISHent.cleanup.build_transcript_gene_map", return_value={"ENS001": "GAPDH"})
def test_annotate_expression_risk_low(
    mock_map, mock_gene, mock_outdir, mock_pcfg, mock_gcfg, task_cleanup
):
    """Low expression (1-10 TPM) should get low risk."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("gene_id,expression\n")
        f.write("ENS001,5.0\n")
        count_table = f.name

    try:
        mock_gcfg.return_value.reference_annotation = "/fake/annotation.gtf"
        mock_pcfg.return_value.encode_count_table = count_table

        df = pd.DataFrame({
            "name": ["p1"],
            "off_target_genes": ["GAPDH(2)"],
        })
        task_cleanup._annotate_expression_risk(df)
        assert "low" in df["expression_risk"].iloc[0]
    finally:
        os.unlink(count_table)


@patch("eFISHent.cleanup.GeneralConfig")
@patch("eFISHent.cleanup.ProbeConfig")
def test_annotate_expression_risk_no_off_target_column(mock_pcfg, mock_gcfg, task_cleanup):
    """Missing off_target_genes column should still work (early return)."""
    mock_gcfg.return_value.reference_annotation = "/fake/annotation.gtf"
    df = pd.DataFrame({"name": ["p1"]})
    task_cleanup._annotate_expression_risk(df)
    assert "expression_risk" in df.columns
    assert df["expression_risk"].iloc[0] == ""


# ---------------------------------------------------------------------------
# _compute_summary
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
def test_compute_summary_basic(mock_outdir, mock_gene, task_cleanup):
    """Basic summary computation with a probe DataFrame."""
    with patch.object(task_cleanup, "_get_gene_length", return_value=5000):
        with patch("eFISHent.console.get_funnel_data", return_value=[("initial", 1000)]):
            df = pd.DataFrame({
                "name": ["p1", "p2", "p3"],
                "length": [20, 22, 24],
                "start": [10, 100, 200],
                "end": [30, 122, 224],
                "TM": [50.0, 52.0, 48.0],
                "GC": [45.0, 50.0, 55.0],
            })
            summary = task_cleanup._compute_summary(df)
            assert summary["probe_count"] == 3
            assert summary["gene_name"] == "ACTB"
            assert summary["initial_count"] == 1000
            assert summary["gene_length"] == 5000
            assert summary["length_range"] == (20, 24)
            assert summary["tm_median"] == pytest.approx(50.0)


@patch("eFISHent.cleanup.util.get_gene_name", return_value="ACTB")
@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
def test_compute_summary_overlapping_intervals(mock_outdir, mock_gene, task_cleanup):
    """Coverage computation should correctly merge overlapping intervals."""
    with patch.object(task_cleanup, "_get_gene_length", return_value=1000):
        with patch("eFISHent.console.get_funnel_data", return_value=[]):
            df = pd.DataFrame({
                "name": ["p1", "p2"],
                "length": [20, 20],
                "start": [10, 20],
                "end": [30, 40],
                "TM": [50.0, 50.0],
                "GC": [48.0, 48.0],
            })
            summary = task_cleanup._compute_summary(df)
            # Overlapping: 10-30 and 20-40 merge to 10-40 = 30 bp
            assert summary["coverage_pct"] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# _add_transcriptome_details
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.GeneralConfig")
def test_add_transcriptome_details_no_txome(mock_gcfg, task_cleanup):
    """Without a reference transcriptome, columns should be initialized to defaults."""
    mock_gcfg.return_value.reference_transcriptome = ""
    mock_gcfg.return_value.reference_annotation = ""
    df = pd.DataFrame({"name": ["p1"], "sequence": ["ATGC"]})
    task_cleanup._add_transcriptome_details(df)
    assert df["txome_off_targets"].iloc[0] == 0
    assert df["off_target_genes"].iloc[0] == ""
    assert df["worst_match"].iloc[0] == ""


# ---------------------------------------------------------------------------
# _blast_final_probes
# ---------------------------------------------------------------------------


@patch("shutil.which", return_value=None)
def test_blast_final_probes_no_blastn(mock_which, task_cleanup):
    """Should return None when blastn is not available."""
    df = pd.DataFrame({"name": ["p1"], "sequence": ["ATGC"]})
    result = task_cleanup._blast_final_probes(df, "/fake/txome.fa")
    assert result is None


# ---------------------------------------------------------------------------
# _run_blast_verification
# ---------------------------------------------------------------------------


@patch("shutil.which", return_value=None)
def test_blast_verification_no_tools(mock_which, task_cleanup):
    """Should return None when BLAST tools are not available."""
    result = task_cleanup._run_blast_verification("/fake/probes.fasta")
    assert result is None


# ---------------------------------------------------------------------------
# _get_gene_length
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.util.get_output_dir", return_value="/fake")
@patch("eFISHent.cleanup.util.get_gene_name", return_value="test")
@patch("eFISHent.cleanup.SequenceConfig")
def test_get_gene_length_no_files(mock_seqcfg, mock_gene, mock_outdir, task_cleanup):
    """Should return 0 when no sequence files exist."""
    mock_seqcfg.return_value.sequence_file = ""
    result = task_cleanup._get_gene_length()
    assert result == 0


@patch("eFISHent.cleanup.util.get_output_dir")
@patch("eFISHent.cleanup.util.get_gene_name", return_value="test")
def test_get_gene_length_from_sequence_file(mock_gene, mock_outdir, task_cleanup):
    """Should read gene length from sequence fasta."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mock_outdir.return_value = tmpdir
        fasta_path = os.path.join(tmpdir, "test_sequence.fasta")
        Bio.SeqIO.write(
            [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("A" * 500), id="test")],
            fasta_path,
            "fasta",
        )
        result = task_cleanup._get_gene_length()
        assert result == 500


# ---------------------------------------------------------------------------
# GC scoring edge cases (internal to _compute_quality_scores)
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_gc_score_optimal_range(mock_cfg):
    """GC in 45-52% should get full GC component score."""
    df_45 = _make_quality_df(GC=[45.0])
    df_52 = _make_quality_df(GC=[52.0])
    q_45 = CleanUpOutput._compute_quality_scores(df_45).iloc[0]
    q_52 = CleanUpOutput._compute_quality_scores(df_52).iloc[0]
    # Both should be roughly equal since GC is in optimal range
    assert abs(q_45 - q_52) < 2.0


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_gc_score_gentle_penalty_range(mock_cfg):
    """GC 52-55% should have gentle penalty, less than >55%."""
    df_53 = _make_quality_df(GC=[53.0])
    df_60 = _make_quality_df(GC=[60.0])
    q_53 = CleanUpOutput._compute_quality_scores(df_53).iloc[0]
    q_60 = CleanUpOutput._compute_quality_scores(df_60).iloc[0]
    assert q_53 > q_60


# ---------------------------------------------------------------------------
# Quality score weight branches (accessibility / binding toggles)
# ---------------------------------------------------------------------------


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_with_accessibility_active(mock_cfg):
    """When accessibility varies (< 1.0), accessibility weight branch is used."""
    df = pd.DataFrame({
        "TM": [50.0, 50.0],
        "GC": [48.0, 48.0],
        "cpg_fraction": [0.02, 0.02],
        "deltaG": [-2.0, -2.0],
        "kmers": [1, 1],
        "count": [0, 0],
        "txome_off_targets": [0, 0],
        "low_complexity": [0.0, 0.0],
        "accessibility": [0.9, 0.3],
        "on_target_dg": [-30.0, -30.0],
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    # Higher accessibility should produce higher score
    assert scores.iloc[0] > scores.iloc[1]


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_with_binding_active(mock_cfg):
    """When on_target_dg varies (some < 0), binding weight branch is used."""
    df = pd.DataFrame({
        "TM": [50.0, 50.0],
        "GC": [48.0, 48.0],
        "cpg_fraction": [0.02, 0.02],
        "deltaG": [-2.0, -2.0],
        "kmers": [1, 1],
        "count": [0, 0],
        "txome_off_targets": [0, 0],
        "low_complexity": [0.0, 0.0],
        "accessibility": [1.0, 1.0],
        "on_target_dg": [-50.0, -10.0],
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    # More negative on_target_dg = more stable binding = higher score
    assert scores.iloc[0] > scores.iloc[1]


@patch("eFISHent.cleanup.ProbeConfig", return_value=_FakeProbeConfig())
def test_quality_score_both_accessibility_and_binding(mock_cfg):
    """When both accessibility and binding are active."""
    df = pd.DataFrame({
        "TM": [50.0, 50.0],
        "GC": [48.0, 48.0],
        "cpg_fraction": [0.02, 0.02],
        "deltaG": [-2.0, -2.0],
        "kmers": [1, 1],
        "count": [0, 0],
        "txome_off_targets": [0, 0],
        "low_complexity": [0.0, 0.0],
        "accessibility": [0.9, 0.3],
        "on_target_dg": [-50.0, -10.0],
    })
    scores = CleanUpOutput._compute_quality_scores(df)
    assert scores.iloc[0] > scores.iloc[1]


# ---------------------------------------------------------------------------
# CSV output formatting (prettify_sequences edge cases)
# ---------------------------------------------------------------------------


def test_prettify_sequences_empty_df(task_cleanup):
    """Empty DataFrame should produce an empty list."""
    df = pd.DataFrame({"name": pd.Series(dtype=str), "sequence": pd.Series(dtype=str)})
    result = task_cleanup.prettify_sequences(df)
    assert result == []


def test_prettify_sequences_preserves_order(task_cleanup):
    """Output should match DataFrame row order."""
    df = pd.DataFrame({
        "name": ["z-probe", "a-probe", "m-probe"],
        "sequence": ["ATGC", "GCGC", "TTTT"],
    })
    result = task_cleanup.prettify_sequences(df)
    assert [r.id for r in result] == ["z-probe", "a-probe", "m-probe"]


# ---------------------------------------------------------------------------
# _GENE_FAMILY_PREFIXES constant validation
# ---------------------------------------------------------------------------


def test_gene_family_prefixes_exist():
    """The gene family prefix list should have known entries."""
    prefixes = CleanUpOutput._GENE_FAMILY_PREFIXES
    assert "ACT" in prefixes
    assert "MYH" in prefixes
    assert "HIST" in prefixes
    assert "HLA" in prefixes
    assert "KRT" in prefixes
    assert "RPL" in prefixes
