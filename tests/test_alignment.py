import os

import pandas as pd
import pytest

from eFISHent.alignment import AlignProbeCandidates
from eFISHent.alignment import PrepareAnnotationFile


def test_prepare_annotations():
    fname_output = "./tests/sacCer3.gtf.parq"
    PrepareAnnotationFile().prepare_gtf_file("./tests/sacCer3.gtf", fname_output)
    assert os.path.isfile(fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "source", "feature", "start", "end"]:
        assert col in df.columns


@pytest.fixture
def task_align() -> AlignProbeCandidates:
    task = AlignProbeCandidates()
    task.fname_genome = "./tests/sacCer3"
    return task


@pytest.mark.parametrize(
    "fname", ["./tests/count_table1.tsv", "./tests/count_table3.txt"]
)
def test_read_count_table(task_align: AlignProbeCandidates, fname):
    task_align.fname_count = fname
    df = task_align.count_table
    for col in ["gene_id", "count"]:
        assert col in df.columns
    assert df["count"].dtype == float


@pytest.mark.parametrize("gene,endogenous", [("renilla", False), ("aad4", True)])
def test_align_probes(task_align: AlignProbeCandidates, gene, endogenous):
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
def test_filter_unique_probes(task_align: AlignProbeCandidates, gene):
    fname_sam = f"./tests/{gene}_basic.sam"
    task_align.fname_gene = gene
    task_align.fname_sam = fname_sam

    df_endo = task_align.filter_unique_probes(is_endogenous=True)
    df_exo = task_align.filter_unique_probes(is_endogenous=False)
    assert "qname" in df_endo
    assert "qname" in df_exo
    assert len(pd.merge(df_endo, df_exo, how="inner", on="qname")) == 0


def test_exclude_gene_of_interest(task_align: AlignProbeCandidates):
    df = pd.read_csv("./tests/count_table2.csv")

    # Using ensembl ID
    ensembl = "ENSG00000281100"
    out = task_align.exclude_gene_of_interest(
        df, ensembl_id=ensembl, fname_full_gene="", threads=2
    )
    assert ensembl not in out["gene_id"]

    # Using blast
    aad4_gene_id = "YDL243C"
    out = task_align.exclude_gene_of_interest(
        df, ensembl_id="", fname_full_gene="./tests/aad4.fasta", threads=2
    )
    assert aad4_gene_id not in out["gene_id"]


def test_join_alignment_with_annotation(task_align: AlignProbeCandidates):
    task_align.fname_sam = "./tests/aad4_basic.sam"
    task_align.fname_count = "./tests/count_table1.tsv"
    task_align.fname_gtf = "./tests/sacCer3.gtf.parq"
    df_sam = task_align.filter_unique_probes(is_endogenous=True)
    df_norm = task_align.norm_table

    df = task_align.join_alignment_with_annotation(df_sam, df_norm)
    assert len(df) >= len(df_sam)
    for col in ["gene_id", "count", "qname"]:
        assert col in df.columns


def test_get_most_expressed_genes(task_align: AlignProbeCandidates):
    task_align.fname_count = "./tests/count_table1.tsv"
    task_align.fname_gtf = "./tests/sacCer3.gtf.parq"

    percentage10 = task_align.get_most_expressed_genes(task_align.norm_table, 10)
    percentage90 = task_align.get_most_expressed_genes(task_align.norm_table, 90)
    assert len(percentage10) >= len(percentage90)
