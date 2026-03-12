import os
import shutil
import tempfile

import pandas as pd
import pytest

from eFISHent.alignment import AlignProbeCandidates

BOWTIE_AVAILABLE = shutil.which("bowtie") is not None
BOWTIE2_AVAILABLE = shutil.which("bowtie2") is not None


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
    df_endo, _ = task_align.filter_unique_probes()
    task_align.is_endogenous = False
    df_exo, _ = task_align.filter_unique_probes()
    assert "qname" in df_endo
    assert "qname" in df_exo
    assert len(pd.merge(df_endo, df_exo, how="inner", on="qname")) == 0


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


@pytest.mark.skipif(not BOWTIE2_AVAILABLE, reason="bowtie2 not installed")
@pytest.mark.parametrize("gene,endogenous", [("renilla", False), ("aad4", True)])
def test_align_probes_bowtie2(gene, endogenous):
    """Verify bowtie2 produces valid SAM output."""
    # Check if bt2 index exists
    if not os.path.isfile("./tests/data/sacCer3.1.bt2"):
        pytest.skip("bowtie2 index not built")

    task = AlignProbeCandidates()
    task.fname_genome = "./tests/data/sacCer3"
    task.aligner = "bowtie2"
    fname_sam = f"./tests/data/{gene}_basic_bt2.sam"
    task.fname_fasta = f"./tests/data/{gene}_basic.fasta"
    task.fname_sam = fname_sam
    task.max_off_targets = 0
    task.is_endogenous = endogenous
    task.align_probes_bowtie2(threads=2)
    assert os.path.isfile(fname_sam)

    # Clean up
    if os.path.exists(fname_sam):
        os.remove(fname_sam)


@pytest.mark.skipif(not BOWTIE2_AVAILABLE, reason="bowtie2 not installed")
def test_filter_unique_probes_bowtie2():
    """Verify filtering works with bowtie2 SAM output."""
    if not os.path.isfile("./tests/data/sacCer3.1.bt2"):
        pytest.skip("bowtie2 index not built")

    task = AlignProbeCandidates()
    task.fname_genome = "./tests/data/sacCer3"
    task.aligner = "bowtie2"
    fname_sam = "./tests/data/aad4_basic_bt2_test.sam"
    task.fname_fasta = "./tests/data/aad4_basic.fasta"
    task.fname_sam = fname_sam
    task.max_off_targets = 0
    task.no_alternative_loci = False

    task.is_endogenous = True
    task.align_probes_bowtie2(threads=2)
    df_endo, _ = task.filter_unique_probes()
    assert "qname" in df_endo.columns

    task.is_endogenous = False
    task.align_probes_bowtie2(threads=2)
    df_exo, _ = task.filter_unique_probes()
    assert "qname" in df_exo.columns

    # Clean up
    if os.path.exists(fname_sam):
        os.remove(fname_sam)


class TestMaskedIntervals:
    """Tests for dustmasker interval parsing and lookup."""

    def test_load_masked_intervals(self):
        """Parse dustmasker interval format correctly."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".intervals", delete=False) as f:
            f.write(">chrI\n0 - 100\n500 - 600\n>chrII\n200 - 300\n")
            f.flush()
            path = f.name

        result = AlignProbeCandidates._load_masked_intervals(path)
        os.unlink(path)

        assert "chrI" in result
        assert "chrII" in result
        assert result["chrI"] == ([0, 500], [100, 600])
        assert result["chrII"] == ([200], [300])

    def test_is_in_masked_region(self):
        """Check masked region overlap detection."""
        task = AlignProbeCandidates()
        task._masked_intervals = {
            "chrI": ([100, 500], [200, 600]),
        }
        # Inside masked region
        assert task._is_in_masked_region("chrI", 150, 20) is True
        # Overlapping boundary
        assert task._is_in_masked_region("chrI", 190, 20) is True
        # Outside masked region
        assert task._is_in_masked_region("chrI", 300, 20) is False
        # Unknown chromosome
        assert task._is_in_masked_region("chrIII", 150, 20) is False
