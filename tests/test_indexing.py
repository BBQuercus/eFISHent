import os

import pandas as pd

from eFISHent.alignment import PrepareAnnotationFile


def test_prepare_annotations():
    fname_output = "./tests/mylittle.gtf.parq"
    PrepareAnnotationFile().prepare_gtf_file("./tests/sacCer3.gtf", fname_output)
    assert os.path.isfile(fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "source", "feature", "start", "end"]:
        assert col in df.columns
