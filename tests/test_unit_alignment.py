import os
import shutil

import pandas as pd
import pytest

from eFISHent.alignment import AlignProbeCandidates

BOWTIE_AVAILABLE = shutil.which("bowtie") is not None


@pytest.fixture
def task_align() -> AlignProbeCandidates:
    task = AlignProbeCandidates()
    task.fname_genome = "./tests/data/sacCer3"
    return task


@pytest.mark.parametrize(
    "fname", ["./tests/data/count_table1.tsv", "./tests/data/count_table3.txt"]
)
def test_read_count_table(task_align: AlignProbeCandidates, fname):
    task_align.fname_count = fname
    df = task_align.count_table
    for col in ["gene_id", "count"]:
        assert col in df.columns
    assert df["count"].dtype == float


@pytest.mark.skipif(not BOWTIE_AVAILABLE, reason="bowtie not installed")
@pytest.mark.parametrize("gene,endogenous", [("renilla", False), ("aad4", True)])
def test_align_probes(task_align: AlignProbeCandidates, gene, endogenous):
    fname_sam = f"./tests/data/{gene}_basic.sam"
    fname_fastq = f"./tests/data/{gene}_basic.fastq"
    task_align.fname_fasta = f"./tests/data/{gene}_basic.fasta"
    task_align.fname_gene = gene
    task_align.fname_sam = fname_sam

    def get_number_of_alignments(fname):
        with open(fname, "r") as f:
            lines = f.readlines()
        lines = list(filter(lambda x: x.startswith("candidate"), lines))
        return len(lines)

    task_align.max_off_targets = 0
    task_align.is_endogenous = endogenous
    task_align.align_probes(threads=2)
    assert os.path.isfile(fname_fastq)
    assert os.path.isfile(fname_sam)

    length1 = get_number_of_alignments(fname_sam)
    task_align.max_off_targets = 20
    task_align.is_endogenous = endogenous
    task_align.align_probes(threads=2)
    length2 = get_number_of_alignments(fname_sam)
    assert length1 <= length2


@pytest.mark.parametrize("gene", ["renilla", "aad4"])
def test_filter_unique_probes(task_align: AlignProbeCandidates, gene):
    fname_sam = f"./tests/data/{gene}_basic.sam"
    task_align.fname_gene = gene
    task_align.fname_sam = fname_sam

    task_align.is_endogenous = True
    df_endo = task_align.filter_unique_probes()
    task_align.is_endogenous = False
    df_exo = task_align.filter_unique_probes()
    assert "qname" in df_endo
    assert "qname" in df_exo
    assert len(pd.merge(df_endo, df_exo, how="inner", on="qname")) == 0


# def test_exclude_gene_of_interest_ensembl(task_align: AlignProbeCandidates):
#     df = pd.read_csv("./tests/data/count_table2.csv")
#     ensembl = "ENSG00000281100"
#     out = task_align.exclude_gene_of_interest(
#         df, ensembl_id=ensembl, fname_full_gene="", threads=2
#     )
#     assert ensembl not in out["gene_id"]


# def test_exclude_gene_of_interest_blast(task_align: AlignProbeCandidates):
#     df = pd.read_csv("./tests/data/count_table2.csv")
#     aad4_gene_id = "YDL243C"
#     out = task_align.exclude_gene_of_interest(
#         df, ensembl_id="", fname_full_gene="./tests/data/aad4.fasta", threads=2
#     )
#     assert aad4_gene_id not in out["gene_id"]


# TODO add test for filter_using_count
# TODO add test for fwd and rev strand
# def test_join_alignment_with_annotation(task_align: AlignProbeCandidates):
#     task_align.fname_sam = "./tests/data/aad4_basic.sam"
#     task_align.fname_count = "./tests/data/count_table1.tsv"
#     task_align.fname_gtf = "./tests/data/sacCer3.gtf.parq"
#     task_align.is_endogenous = True
#     df_sam = task_align.filter_unique_probes()
#     df_norm = task_align.norm_table

#     df = task_align.join_alignment_with_annotation(df_sam, df_norm)
#     assert len(df) >= len(df_sam)
#     for col in ["gene_id", "count", "qname"]:
#         assert col in df.columns


# def test_get_most_expressed_genes(task_align: AlignProbeCandidates):
#     task_align.fname_count = "./tests/data/count_table1.tsv"
#     task_align.fname_gtf = "./tests/data/sacCer3.gtf.parq"

#     percentage10 = task_align.get_most_expressed_genes(task_align.norm_table, 10)
#     percentage90 = task_align.get_most_expressed_genes(task_align.norm_table, 90)
#     assert len(percentage10) >= len(percentage90)


class TestParseRawPysam:
    """Tests for parse_raw_pysam static method."""

    def test_empty_input(self):
        """Empty input should return empty DataFrame with correct columns."""
        result = AlignProbeCandidates.parse_raw_pysam("")
        assert len(result) == 0
        expected_cols = [
            "qname", "flag", "rname", "pos", "mapq",
            "cigar", "rnext", "pnext", "tlen", "seq"
        ]
        assert list(result.columns) == expected_cols

    def test_single_alignment(self):
        """Single alignment row should parse correctly."""
        # SAM format: qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, ...
        sam_line = "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        assert len(result) == 1
        assert result.iloc[0]["qname"] == "candidate-0-100"
        assert result.iloc[0]["flag"] == "0"
        assert result.iloc[0]["rname"] == "chrI"
        assert result.iloc[0]["pos"] == "1000"
        assert result.iloc[0]["mapq"] == "60"
        assert result.iloc[0]["cigar"] == "21M"

    def test_multiple_alignments(self):
        """Multiple alignment rows should all be parsed."""
        sam_lines = "\n".join([
            "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*",
            "candidate-1-150\t0\tchrII\t2000\t60\t22M\t*\t0\t0\tATGCATGCATGCATGCATGCAT\t*",
            "candidate-2-200\t16\tchrIII\t3000\t60\t23M\t*\t0\t0\tATGCATGCATGCATGCATGCATG\t*",
        ])
        result = AlignProbeCandidates.parse_raw_pysam(sam_lines)

        assert len(result) == 3
        assert list(result["qname"]) == ["candidate-0-100", "candidate-1-150", "candidate-2-200"]
        assert list(result["rname"]) == ["chrI", "chrII", "chrIII"]

    def test_extra_columns_ignored(self):
        """Extra columns beyond the first 10 should be ignored."""
        # SAM can have optional fields (XA, MD, NM, etc.)
        sam_line = "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGC\t*\tXA:Z:extra\tMD:Z:21\tNM:i:0"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        assert len(result) == 1
        assert len(result.columns) == 10
        assert "XA" not in result.columns

    def test_trailing_newline_handled(self):
        """Trailing newlines should not create empty rows."""
        sam_line = "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGC\t*\n"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        # dropna() should remove the incomplete row from trailing newline
        assert len(result) == 1

    def test_unmapped_read(self):
        """Unmapped reads (flag 4) should be parsed correctly."""
        sam_line = "candidate-0-100\t4\t*\t0\t0\t*\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        assert len(result) == 1
        assert result.iloc[0]["flag"] == "4"
        assert result.iloc[0]["rname"] == "*"

    def test_reverse_strand(self):
        """Reverse strand alignments (flag 16) should be parsed correctly."""
        sam_line = "candidate-0-100\t16\tchrI\t1000\t60\t21M\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        assert len(result) == 1
        assert result.iloc[0]["flag"] == "16"
