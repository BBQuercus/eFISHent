import os

import pytest
import pandas as pd

from eFISHent.indexing import PrepareAnnotationFile


@pytest.fixture
def task_prepare():
    return PrepareAnnotationFile()


def test_prepare_annotations(task_prepare):
    fname_output = "./mylittle.gtf.parq"
    task_prepare.prepare_gtf_file("./tests/data/sacCer3.gtf", fname_output)

    assert os.path.isfile(fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "source", "feature", "start", "end"]:
        assert col in df.columns
    os.remove(fname_output)
