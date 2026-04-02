import shutil
import tempfile

import pandas as pd
import luigi
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import pytest

from eFISHent.transcriptome_filter import TranscriptomeFiltering

BLAST_AVAILABLE = shutil.which("blastn") is not None


def test_blast_output_parsing():
    """Parse sample BLAST -outfmt 6 output."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "candidate-0-100\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "candidate-0-100\tGENE2\t80.0\t18\t3\t1\t1\t18\t200\t217\t0.01\t25\n"
        "candidate-1-200\tGENE1\t95.0\t22\t1\t0\t1\t22\t300\t321\t0.0001\t35\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)
    assert len(df) == 3
    assert df.iloc[0]["qseqid"] == "candidate-0-100"
    assert df.iloc[0]["pident"] == 90.0


def test_identity_threshold_filtering():
    """Hits below threshold are ignored."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "probe1\tGENE1\t90.0\t20\t2\t0\t1\t20\t100\t119\t0.001\t30\n"
        "probe1\tGENE2\t60.0\t18\t7\t1\t1\t18\t200\t217\t0.5\t10\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)

    # Apply threshold at 75%
    threshold = 75.0
    min_match_len = 15
    df["effective_len"] = df["length"] - df["gapopen"]
    df_hits = df[
        (df["effective_len"] >= min_match_len) & (df["pident"] >= threshold)
    ]
    assert len(df_hits) == 1
    assert df_hits.iloc[0]["sseqid"] == "GENE1"


def test_target_gene_excluded():
    """Self-hits to target gene don't count as off-targets."""
    from io import StringIO

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    sample = (
        "probe1\trenilla_transcript\t100.0\t25\t0\t0\t1\t25\t1\t25\t0.0\t50\n"
        "probe1\tother_gene\t85.0\t20\t3\t0\t1\t20\t100\t119\t0.01\t30\n"
    )
    df = pd.read_csv(StringIO(sample), sep="\t", header=None, names=columns)

    gene_name = "renilla"
    off_targets = df[~df["sseqid"].str.lower().str.contains(gene_name, na=False)]
    assert len(off_targets) == 1
    assert off_targets.iloc[0]["sseqid"] == "other_gene"


def test_cross_hybridization_probe_removed(monkeypatch):
    """A >=16nt, >=95% identity off-target transcript hit should veto the probe."""
    cfg = luigi.configuration.get_config()
    old_values = {}
    config_updates = {
        ("GeneralConfig", "reference_transcriptome"): "/tmp/fake_txome.fa",
        ("GeneralConfig", "reference_annotation"): "",
        ("GeneralConfig", "threads"): "1",
        ("GeneralConfig", "output_dir"): "",
        ("ProbeConfig", "max_transcriptome_off_targets"): "10",
        ("ProbeConfig", "blast_identity_threshold"): "75.0",
        ("ProbeConfig", "min_blast_match_length"): "0",
        ("ProbeConfig", "min_length"): "20",
        ("ProbeConfig", "max_length"): "25",
        ("ProbeConfig", "reject_cross_hybridization"): "true",
    }

    for (section, key), value in config_updates.items():
        old_values[(section, key)] = cfg.get(section, key, fallback=None)
        if not cfg.has_section(section):
            cfg.add_section(section)
        cfg.set(section, key, value)

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            probe_fasta = f"{tmpdir}/probes.fasta"
            output_fasta = f"{tmpdir}/out.fasta"
            output_csv = f"{tmpdir}/hits.csv"
            blast_out = f"{tmpdir}/probes_blast.tsv"

            records = [
                Bio.SeqRecord.SeqRecord(
                    Bio.Seq.Seq("ATGCGTACGTAGCTAGCTAGC"),
                    id="cross-hyb", name="cross-hyb", description=""
                ),
                Bio.SeqRecord.SeqRecord(
                    Bio.Seq.Seq("TTTTGGGGAAAACCCCTTTTG"),
                    id="safe-probe", name="safe-probe", description=""
                ),
            ]
            Bio.SeqIO.write(records, probe_fasta, "fasta")

            task = TranscriptomeFiltering()
            monkeypatch.setattr(
                task,
                "input",
                lambda: {
                    "probes": {"fasta": luigi.LocalTarget(probe_fasta)},
                    "blastdb": luigi.LocalTarget("/tmp/fake_txome.fa.nsq"),
                },
            )
            monkeypatch.setattr(
                task,
                "output",
                lambda: {
                    "fasta": luigi.LocalTarget(output_fasta),
                    "table": luigi.LocalTarget(output_csv),
                },
            )
            monkeypatch.setattr("eFISHent.transcriptome_filter.util.get_gene_name", lambda hashed=False: "targetgene")
            monkeypatch.setattr("eFISHent.transcriptome_filter.util.log_stage_start", lambda *args, **kwargs: None)
            monkeypatch.setattr("eFISHent.transcriptome_filter.util.log_and_check_candidates", lambda *args, **kwargs: None)

            def fake_check_call(args, stdout=None, stderr=None):
                assert args[0] == "blastn"
                with open(blast_out, "w") as handle:
                    handle.write(
                        "cross-hyb\tother_gene\t95.0\t16\t0\t0\t1\t16\t100\t115\t1e-5\t40\n"
                    )
                    handle.write(
                        "safe-probe\ttargetgene_transcript\t100.0\t21\t0\t0\t1\t21\t1\t21\t0.0\t50\n"
                    )
                return 0

            monkeypatch.setattr("eFISHent.transcriptome_filter.subprocess.check_call", fake_check_call)

            task.run()

            kept = [record.id for record in Bio.SeqIO.parse(output_fasta, "fasta")]
            assert kept == ["safe-probe"]
    finally:
        for (section, key), old_value in old_values.items():
            if old_value is None:
                cfg.remove_option(section, key)
            else:
                cfg.set(section, key, old_value)
