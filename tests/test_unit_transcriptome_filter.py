import os
import shutil
import tempfile

import pandas as pd
import luigi
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

from eFISHent.transcriptome_filter import TranscriptomeFiltering, BuildTranscriptomeBlastDB

BLAST_AVAILABLE = shutil.which("blastn") is not None

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
]


def _make_records(*specs):
    """Create SeqRecord objects from (id, seq) tuples."""
    return [
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(seq), id=name, name=name, description=""
        )
        for name, seq in specs
    ]


def _luigi_config_ctx(updates: dict):
    """Context-manager-like helper returning (old_values, cfg) for luigi config.

    Usage:
        old, cfg = _luigi_config_ctx({("Section", "key"): "val"})
        ...  # do work
        _luigi_config_restore(old, cfg)
    """
    cfg = luigi.configuration.get_config()
    old_values = {}
    for (section, key), value in updates.items():
        old_values[(section, key)] = cfg.get(section, key, fallback=None)
        if not cfg.has_section(section):
            cfg.add_section(section)
        cfg.set(section, key, value)
    return old_values, cfg


def _luigi_config_restore(old_values, cfg):
    for (section, key), old_value in old_values.items():
        if old_value is None:
            cfg.remove_option(section, key)
        else:
            cfg.set(section, key, old_value)


# Shared luigi config overrides used by most integration-style tests
_BASE_CONFIG = {
    ("GeneralConfig", "reference_transcriptome"): "/tmp/fake_txome.fa",
    ("GeneralConfig", "reference_annotation"): "",
    ("GeneralConfig", "threads"): "1",
    ("GeneralConfig", "output_dir"): "",
    ("ProbeConfig", "max_transcriptome_off_targets"): "0",
    ("ProbeConfig", "blast_identity_threshold"): "75.0",
    ("ProbeConfig", "min_blast_match_length"): "0",
    ("ProbeConfig", "min_length"): "20",
    ("ProbeConfig", "max_length"): "25",
    ("ProbeConfig", "reject_cross_hybridization"): "true",
}


def _setup_task(monkeypatch, tmpdir, records, blast_lines, *, gene_name="targetgene", extra_config=None):
    """Wire up a TranscriptomeFiltering task with mocked I/O and BLAST."""
    config = dict(_BASE_CONFIG)
    if extra_config:
        config.update(extra_config)

    old, cfg = _luigi_config_ctx(config)

    probe_fasta = os.path.join(tmpdir, "probes.fasta")
    output_fasta = os.path.join(tmpdir, "out.fasta")
    output_csv = os.path.join(tmpdir, "hits.csv")
    blast_out = os.path.join(tmpdir, "probes_blast.tsv")

    Bio.SeqIO.write(records, probe_fasta, "fasta")

    task = TranscriptomeFiltering()
    monkeypatch.setattr(
        task, "input",
        lambda: {
            "probes": {"fasta": luigi.LocalTarget(probe_fasta)},
            "blastdb": luigi.LocalTarget("/tmp/fake_txome.fa.nsq"),
        },
    )
    monkeypatch.setattr(
        task, "output",
        lambda: {
            "fasta": luigi.LocalTarget(output_fasta),
            "table": luigi.LocalTarget(output_csv),
        },
    )
    monkeypatch.setattr(
        "eFISHent.transcriptome_filter.util.get_gene_name",
        lambda hashed=False: gene_name,
    )
    monkeypatch.setattr(
        "eFISHent.transcriptome_filter.util.log_stage_start",
        lambda *a, **kw: None,
    )
    monkeypatch.setattr(
        "eFISHent.transcriptome_filter.util.log_and_check_candidates",
        lambda *a, **kw: None,
    )

    def fake_check_call(args, stdout=None, stderr=None):
        assert args[0] == "blastn"
        with open(blast_out, "w") as fh:
            fh.write(blast_lines)
        return 0

    monkeypatch.setattr(
        "eFISHent.transcriptome_filter.subprocess.check_call", fake_check_call,
    )

    return task, output_fasta, output_csv, old, cfg


# ---------------------------------------------------------------------------
# 1. BLAST output parsing
# ---------------------------------------------------------------------------


def test_blast_output_parsing():
    """Parse sample BLAST -outfmt 6 output."""
    from io import StringIO

    sample = (
        "candidate-0-100\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "candidate-0-100\tGENE2\t80.0\t18\t3\t1\t1\t18\t200\t217\t0.01\t25\n"
        "candidate-1-200\tGENE1\t95.0\t22\t1\t0\t1\t22\t300\t321\t0.0001\t35\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=BLAST_COLUMNS)
    assert len(df) == 3
    assert df.iloc[0]["qseqid"] == "candidate-0-100"
    assert df.iloc[0]["pident"] == 90.0


def test_blast_output_empty_results():
    """Empty BLAST result file produces an empty DataFrame."""
    from io import StringIO

    df = pd.read_csv(StringIO(""), sep="\t", header=None, names=BLAST_COLUMNS)
    assert len(df) == 0
    assert list(df.columns) == BLAST_COLUMNS


def test_blast_output_single_row():
    """Single BLAST hit parses correctly."""
    from io import StringIO

    sample = "probe1\tGENE1\t99.5\t25\t0\t0\t1\t25\t1\t25\t1e-10\t50\n"
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=BLAST_COLUMNS)
    assert len(df) == 1
    assert df.iloc[0]["bitscore"] == 50.0


# ---------------------------------------------------------------------------
# 2. Identity threshold filtering
# ---------------------------------------------------------------------------


def test_identity_threshold_filtering():
    """Hits below threshold are ignored."""
    from io import StringIO

    sample = (
        "probe1\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "probe1\tGENE2\t60.0\t18\t7\t1\t1\t18\t200\t217\t0.5\t10\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=BLAST_COLUMNS)
    threshold = 75.0
    min_match_len = 15
    df["effective_len"] = df["length"] - df["gapopen"]
    df_hits = df[
        (df["effective_len"] >= min_match_len) & (df["pident"] >= threshold)
    ]
    assert len(df_hits) == 1
    assert df_hits.iloc[0]["sseqid"] == "GENE1"


def test_effective_length_filter():
    """Hits with short effective length (alignment - gaps) are excluded."""
    from io import StringIO

    sample = (
        "probe1\tGENE1\t90.0\t20\t2\t3\t1\t20\t100\t119\t0.001\t30\n"  # eff=17
        "probe1\tGENE2\t90.0\t22\t1\t0\t1\t22\t200\t221\t0.001\t30\n"  # eff=22
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=BLAST_COLUMNS)
    df["effective_len"] = df["length"] - df["gapopen"]
    min_match_len = 18
    df_hits = df[
        (df["effective_len"] >= min_match_len) & (df["pident"] >= 75.0)
    ]
    assert len(df_hits) == 1
    assert df_hits.iloc[0]["sseqid"] == "GENE2"


# ---------------------------------------------------------------------------
# 3. Self-hit exclusion
# ---------------------------------------------------------------------------


def test_target_gene_excluded():
    """Self-hits to target gene don't count as off-targets."""
    from io import StringIO

    sample = (
        "probe1\trenilla_transcript\t100.0\t25\t0\t0\t1\t25\t1\t25\t0.0\t50\n"
        "probe1\tother_gene\t85.0\t20\t3\t0\t1\t20\t100\t119\t0.01\t30\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=BLAST_COLUMNS)
    gene_name = "renilla"
    off_targets = df[~df["sseqid"].str.lower().str.contains(gene_name, na=False)]
    assert len(off_targets) == 1
    assert off_targets.iloc[0]["sseqid"] == "other_gene"


# ---------------------------------------------------------------------------
# 4. Full run() tests via mocked subprocess
# ---------------------------------------------------------------------------


def test_probe_no_blast_hits_passes(monkeypatch):
    """Probe with no BLAST hits passes filter."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=""
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


def test_probe_self_hit_passes(monkeypatch):
    """Probe with hit only to self-gene passes filter."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    blast = "probe1\ttargetgene_transcript\t100.0\t21\t0\t0\t1\t21\t1\t21\t0.0\t50\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


def test_probe_off_target_above_threshold_filtered(monkeypatch):
    """Probe with off-target hit above threshold is filtered out."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Hit to other gene: 90% identity, length 20, no gaps → passes both filters
    blast = "probe1\tother_gene\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "0"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


def test_probe_hit_below_identity_passes(monkeypatch):
    """Probe with hit below identity threshold passes filter."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # 60% identity is below threshold of 75%
    blast = "probe1\tother_gene\t60.0\t20\t8\t0\t1\t20\t100\t119\t1.0\t10\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


def test_cross_hybridization_probe_removed(monkeypatch):
    """A >=16nt, >=95% identity off-target transcript hit should veto the probe."""
    records = _make_records(
        ("cross-hyb", "ATGCGTACGTAGCTAGCTAGC"),
        ("safe-probe", "TTTTGGGGAAAACCCCTTTTG"),
    )
    blast = (
        "cross-hyb\tother_gene\t95.0\t16\t0\t0\t1\t16\t100\t115\t1e-5\t40\n"
        "safe-probe\ttargetgene_transcript\t100.0\t21\t0\t0\t1\t21\t1\t21\t0.0\t50\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["safe-probe"]
        finally:
            _luigi_config_restore(old, cfg)


def test_cross_hybridization_warning_only(monkeypatch):
    """With reject_cross_hybridization=false, probes are warned but kept."""
    records = _make_records(
        ("cross-hyb", "ATGCGTACGTAGCTAGCTAGC"),
    )
    blast = "cross-hyb\tother_gene\t96.0\t18\t0\t0\t1\t18\t100\t117\t1e-6\t45\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={
                ("ProbeConfig", "reject_cross_hybridization"): "false",
                ("ProbeConfig", "max_transcriptome_off_targets"): "10",
            },
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # Probe kept because we only warn, not reject
            assert kept == ["cross-hyb"]
        finally:
            _luigi_config_restore(old, cfg)


def test_empty_probe_set(monkeypatch):
    """Empty probe set produces empty output."""
    records = []
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines="",
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


def test_all_probes_filtered(monkeypatch):
    """All probes filtered out → empty output."""
    records = _make_records(
        ("probe1", "ATGCGTACGTAGCTAGCTAGC"),
        ("probe2", "CCCCAAAATTTTGGGGCCCCA"),
    )
    blast = (
        "probe1\toff_target_A\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "probe2\toff_target_B\t88.0\t21\t3\t0\t1\t21\t200\t220\t0.002\t28\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "0"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


def test_single_probe_passes(monkeypatch):
    """Single probe with no off-target passes."""
    records = _make_records(("only-probe", "ATGCGTACGTAGCTAGCTAGC"))
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines="",
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["only-probe"]
        finally:
            _luigi_config_restore(old, cfg)


def test_multiple_off_targets_counted(monkeypatch):
    """Multiple distinct off-targets are counted independently per probe."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Two distinct off-target genes
    blast = (
        "probe1\toff_A\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "probe1\toff_B\t85.0\t21\t3\t0\t1\t21\t200\t220\t0.01\t25\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        # Allow 1 off-target, but probe has 2 → filtered
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "1"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


def test_off_targets_within_limit_passes(monkeypatch):
    """Probe with off-targets within the allowed limit passes."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    blast = "probe1\toff_A\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        # Allow 1 off-target, probe has exactly 1 → passes
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "1"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


# ---------------------------------------------------------------------------
# 5. Gene name cleaning / suffix stripping
# ---------------------------------------------------------------------------


def test_cleaned_gene_name_self_hit(monkeypatch):
    """Gene name with suffix (e.g. 'EIF2B1_cds') is stripped for self-hit matching."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Subject ID contains cleaned gene name "eif2b1" (gene_name = "eif2b1_cds")
    blast = "probe1\tENST00000|eif2b1_variant\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            gene_name="eif2b1_cds",
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # Should pass because eif2b1 (cleaned from eif2b1_cds) matches
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


# ---------------------------------------------------------------------------
# 6. GTF-based transcript-to-gene mapping
# ---------------------------------------------------------------------------


def test_gtf_transcript_mapping_self_hit(monkeypatch):
    """Self-hit detection via GTF transcript-to-gene mapping."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Use an Ensembl-style subject ID not containing the gene name
    blast = "probe1\tENST00000123456.5\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        # Provide a fake GTF path to trigger the mapping branch
        fake_gtf = os.path.join(tmpdir, "fake.gtf")
        with open(fake_gtf, "w") as fh:
            fh.write("")  # content doesn't matter, we mock the function

        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            gene_name="mygene",
            extra_config={("GeneralConfig", "reference_annotation"): fake_gtf},
        )
        # Mock build_transcript_gene_map at its source module
        monkeypatch.setattr(
            "eFISHent.gene_annotation.build_transcript_gene_map",
            lambda path: {"ENST00000123456": "mygene", "ENST00000999999": "othergene"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # ENST00000123456.5 → base ENST00000123456 → maps to mygene → self-hit
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


def test_gtf_mapping_off_target_detected(monkeypatch):
    """GTF mapping: hit to transcript of different gene is off-target."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    blast = "probe1\tENST00000999999\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        fake_gtf = os.path.join(tmpdir, "fake.gtf")
        with open(fake_gtf, "w") as fh:
            fh.write("")

        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            gene_name="mygene",
            extra_config={
                ("GeneralConfig", "reference_annotation"): fake_gtf,
                ("ProbeConfig", "max_transcriptome_off_targets"): "0",
            },
        )
        monkeypatch.setattr(
            "eFISHent.gene_annotation.build_transcript_gene_map",
            lambda path: {"ENST00000123456": "mygene", "ENST00000999999": "othergene"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


def test_gtf_mapping_failure_falls_back(monkeypatch):
    """When GTF mapping raises, fall back to name-based exclusion."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Subject contains gene name → name-based self-hit detection works
    blast = "probe1\tmygene_transcript\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        fake_gtf = os.path.join(tmpdir, "fake.gtf")
        with open(fake_gtf, "w") as fh:
            fh.write("")

        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            gene_name="mygene",
            extra_config={("GeneralConfig", "reference_annotation"): fake_gtf},
        )

        def raise_mapping(path):
            raise RuntimeError("GTF parse error")

        monkeypatch.setattr(
            "eFISHent.gene_annotation.build_transcript_gene_map",
            raise_mapping,
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # Falls back to name matching → self-hit → probe passes
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


# ---------------------------------------------------------------------------
# 7. BuildTranscriptomeBlastDB
# ---------------------------------------------------------------------------


def test_build_blast_db_output(monkeypatch):
    """BuildTranscriptomeBlastDB.output() returns path with .nsq suffix."""
    cfg = luigi.configuration.get_config()
    old_val = cfg.get("GeneralConfig", "reference_transcriptome", fallback=None)
    if not cfg.has_section("GeneralConfig"):
        cfg.add_section("GeneralConfig")
    cfg.set("GeneralConfig", "reference_transcriptome", "/tmp/my_txome.fa")
    try:
        task = BuildTranscriptomeBlastDB()
        out = task.output()
        assert out.path.endswith(".nsq")
        assert "my_txome.fa" in out.path
    finally:
        if old_val is None:
            cfg.remove_option("GeneralConfig", "reference_transcriptome")
        else:
            cfg.set("GeneralConfig", "reference_transcriptome", old_val)


def test_build_blast_db_run(monkeypatch):
    """BuildTranscriptomeBlastDB.run() calls makeblastdb."""
    cfg = luigi.configuration.get_config()
    old_val = cfg.get("GeneralConfig", "reference_transcriptome", fallback=None)
    if not cfg.has_section("GeneralConfig"):
        cfg.add_section("GeneralConfig")
    cfg.set("GeneralConfig", "reference_transcriptome", "/tmp/my_txome.fa")

    called_with = {}

    def fake_check_call(args, stdout=None, stderr=None):
        called_with["args"] = args
        return 0

    monkeypatch.setattr(
        "eFISHent.transcriptome_filter.subprocess.check_call", fake_check_call,
    )
    try:
        task = BuildTranscriptomeBlastDB()
        task.run()
        assert called_with["args"][0] == "makeblastdb"
        assert "-in" in called_with["args"]
    finally:
        if old_val is None:
            cfg.remove_option("GeneralConfig", "reference_transcriptome")
        else:
            cfg.set("GeneralConfig", "reference_transcriptome", old_val)


# ---------------------------------------------------------------------------
# 8. Min blast match length auto-calculation
# ---------------------------------------------------------------------------


def test_min_blast_match_length_auto(monkeypatch):
    """When min_blast_match_length=0, it auto-calculates from min_length."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Hit with effective length 15 (below auto-threshold of max(18, 0.8*20)=18)
    blast = "probe1\toff_target\t90.0\t15\t1\t0\t1\t15\t100\t114\t0.01\t20\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "min_blast_match_length"): "0"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # Hit is 15nt, auto threshold is 18 → hit filtered → probe passes
            assert kept == ["probe1"]
        finally:
            _luigi_config_restore(old, cfg)


def test_min_blast_match_length_explicit(monkeypatch):
    """Explicit min_blast_match_length is used when set."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    # Hit with effective length 15, explicit threshold 14 → hit counts
    blast = "probe1\toff_target\t90.0\t15\t1\t0\t1\t15\t100\t114\t0.01\t20\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={
                ("ProbeConfig", "min_blast_match_length"): "14",
                ("ProbeConfig", "max_transcriptome_off_targets"): "0",
            },
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            # eff_len=15 >= 14, pident=90 >= 75 → off-target → filtered
            assert kept == []
        finally:
            _luigi_config_restore(old, cfg)


# ---------------------------------------------------------------------------
# 9. Hits table output
# ---------------------------------------------------------------------------


def test_hits_table_saved(monkeypatch):
    """The hits CSV is written with the right columns."""
    records = _make_records(("probe1", "ATGCGTACGTAGCTAGCTAGC"))
    blast = "probe1\toff_target\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
    with tempfile.TemporaryDirectory() as tmpdir:
        task, out_fa, out_csv, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "10"},
        )
        try:
            task.run()
            df = pd.read_csv(out_csv)
            assert "qseqid" in df.columns
            assert "effective_len" in df.columns
            assert len(df) == 1
        finally:
            _luigi_config_restore(old, cfg)


# ---------------------------------------------------------------------------
# 10. Mixed probes — some pass, some fail
# ---------------------------------------------------------------------------


def test_mixed_probes_partial_filtering(monkeypatch):
    """Only probes exceeding off-target limit are filtered."""
    records = _make_records(
        ("good", "ATGCGTACGTAGCTAGCTAGC"),
        ("bad", "CCCCAAAATTTTGGGGCCCCA"),
        ("ok", "GGGGTTTTAAAACCCCGGGGT"),
    )
    blast = (
        "good\ttargetgene_tx\t100.0\t21\t0\t0\t1\t21\t1\t21\t0.0\t50\n"
        "bad\toff_A\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "bad\toff_B\t88.0\t21\t3\t0\t1\t21\t200\t220\t0.002\t28\n"
        "ok\toff_C\t85.0\t20\t3\t0\t1\t20\t300\t319\t0.01\t25\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        # Allow 1 off-target: good(0)→pass, bad(2)→fail, ok(1)→pass
        task, out_fa, _, old, cfg = _setup_task(
            monkeypatch, tmpdir, records, blast_lines=blast,
            extra_config={("ProbeConfig", "max_transcriptome_off_targets"): "1"},
        )
        try:
            task.run()
            kept = [r.id for r in Bio.SeqIO.parse(out_fa, "fasta")]
            assert kept == ["good", "ok"]
        finally:
            _luigi_config_restore(old, cfg)
