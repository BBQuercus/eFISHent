import os
import shutil
import tempfile

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
