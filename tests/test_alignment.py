import os

import pandas as pd
import pytest

from eFISHent.alignment import AlignProbeCandidates
from eFISHent.alignment import PrepareAnnotationFile
from eFISHent.constants import COUNTS_COLUMNS


def test_prepare_annotations():
    fname_output = "./tests/sacCer3.gtf.parq"
    PrepareAnnotationFile().prepare_gtf_file("./tests/sacCer3.gtf", fname_output)
    assert os.path.isfile(fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "source", "feature", "start", "end"]:
        assert col in df.columns

    os.remove(fname_output)


@pytest.fixture
def task_align():
    task = AlignProbeCandidates()
    task.fname_genome = "./tests/sacCer3"
    return task


def test_read_count_table(task_align):
    fname_counts = "./tests/count_table1.tsv"
    df = task_align.read_count_table(fname_counts)
    for col in COUNTS_COLUMNS:
        assert col in df.columns

    # Bad count table with clean_id instead of gene_id
    fname_bad = "./tests/count_table2.csv"
    with pytest.raises(ValueError):
        task_align.read_count_table(fname_bad)


@pytest.mark.parametrize("gene,endogenous", [("renilla", False), ("aad4", True)])
def test_align_probes(task_align, gene, endogenous):
    fname_sam = f"./tests/{gene}_basic.sam"
    fname_fastq = f"./tests/{gene}_basic.fastq"
    task_align.fname_fasta = f"./tests/{gene}_basic.fasta"
    task_align.fname_gene = gene
    task_align.fname_sam = fname_sam

    def get_number_of_alignments(fname):
        with open(fname, "r") as f:
            lines = f.readlines()
        lines = list(filter(lambda x: x.startswith("candidate"), lines))
        return len(lines)

    task_align.align_probes(max_off_targets=0, is_endogenous=endogenous, threads=2)
    assert os.path.isfile(fname_fastq)
    assert os.path.isfile(fname_sam)

    length1 = get_number_of_alignments(fname_sam)
    task_align.align_probes(max_off_targets=20, is_endogenous=endogenous, threads=2)
    length2 = get_number_of_alignments(fname_sam)
    assert length1 <= length2


@pytest.mark.parametrize("gene", ["renilla", "aad4"])
def test_filter_unique_probes(task_align, gene):
    fname_sam = f"./tests/{gene}_basic.sam"
    task_align.fname_gene = gene
    task_align.fname_sam = fname_sam

    df_endo = task_align.filter_unique_probes(is_endogenous=True)
    df_exo = task_align.filter_unique_probes(is_endogenous=False)
    assert "qname" in df_endo
    assert "qname" in df_exo
    assert len(pd.merge(df_endo, df_exo, how="inner", on="qname")) == 0


def test_exclude_gene_of_interest(task_align):
    df = pd.read_csv("./tests/count_table2.csv")

    # Using ensembl ID
    ensembl = "ENSG00000281100"
    df = task_align.exclude_gene_of_interest(
        df, ensembl_id=ensembl, fname_full_gene="", threads=2
    )
    assert ensembl not in df["clean_id"]

    # Using blast
    # Create count_table with gtf merge!
    # aad4_gene_id = "YDL243C"
    # df = task_align.filter_gene_of_interest(
    #     df, ensembl_id="", fname_full_gene="./tests/aad4.fasta", threads=2
    # )
    # assert aad4_gene_id not in df["clean_id"]


def test_get_maximum_fpkm(task_align):
    # task_align.get_maximum_fpkm(df_sam, df_counts, df_gtf)
    pass
