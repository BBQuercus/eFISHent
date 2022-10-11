import Bio.Seq
import Bio.SeqRecord
import luigi
import pandas as pd
import pytest

from eFISHent.cleanup import CleanUpOutput


@pytest.fixture
def task_cleanup():
    return CleanUpOutput()


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
    output = task_cleanup.prettify_table(
        sequences,
        jellyfish_path="./tests/data/sacCer3_15.jf",
        basename="test",
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
    ]
    for col in columns:
        assert col in output.columns
    assert len(output) == len(sequences)
    assert output.loc[0, "name"] == "test-1"
    assert output.loc[len(sequences) - 1, "name"] == f"test-{len(sequences)}"
    assert output.loc[0, "GC"] == 50
    assert output.loc[3, "GC"] == 66.67


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
