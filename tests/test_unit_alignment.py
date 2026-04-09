import os
import shutil
import tempfile
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from eFISHent.alignment import AlignProbeCandidates

BOWTIE_AVAILABLE = shutil.which("bowtie") is not None
BOWTIE2_AVAILABLE = shutil.which("bowtie2") is not None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_sam_line(qname, flag, rname, pos, mapq, seq="ATGCATGCATGCATGCATGCA"):
    """Build a single SAM-format line."""
    cigar = f"{len(seq)}M"
    return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*"


def _make_task(**overrides):
    """Create an AlignProbeCandidates with sensible defaults, no Luigi config needed."""
    task = AlignProbeCandidates()
    defaults = dict(
        fname_genome="./tests/data/sacCer3",
        is_endogenous=True,
        max_off_targets=0,
        no_alternative_loci=False,
        aligner="bowtie2",
        encode_count_table="",
        off_target_min_tm=0,
        mask_repeats=False,
        intergenic_off_targets=False,
        filter_rrna=False,
        has_expression_filter=False,
        has_annotation=False,
        has_intergenic_filter=False,
        has_rrna_filter=False,
        has_pseudogene_filter=False,
        has_transcriptome=False,
    )
    defaults.update(overrides)
    for k, v in defaults.items():
        setattr(task, k, v)
    return task


# ---------------------------------------------------------------------------
# Fixture (original)
# ---------------------------------------------------------------------------

@pytest.fixture
def task_align() -> AlignProbeCandidates:
    task = AlignProbeCandidates()
    task.fname_genome = "./tests/data/sacCer3"
    return task


# ===========================================================================
# 1. Count table parsing (original tests preserved)
# ===========================================================================

@pytest.mark.parametrize(
    "fname", ["./tests/data/count_table1.tsv", "./tests/data/count_table3.txt"]
)
def test_read_count_table(task_align: AlignProbeCandidates, fname):
    task_align.fname_count = fname
    df = task_align.count_table
    for col in ["gene_id", "count"]:
        assert col in df.columns
    assert df["count"].dtype == float


# ===========================================================================
# 2. Original integration tests (require tools installed)
# ===========================================================================

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


# ===========================================================================
# 3. parse_raw_pysam (original tests preserved)
# ===========================================================================

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
        sam_line = "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGC\t*\tXA:Z:extra\tMD:Z:21\tNM:i:0"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

        assert len(result) == 1
        assert len(result.columns) == 10
        assert "XA" not in result.columns

    def test_trailing_newline_handled(self):
        """Trailing newlines should not create empty rows."""
        sam_line = "candidate-0-100\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGC\t*\n"
        result = AlignProbeCandidates.parse_raw_pysam(sam_line)

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


# ===========================================================================
# 4. Original bowtie2 integration tests (preserved)
# ===========================================================================

@pytest.mark.skipif(not BOWTIE2_AVAILABLE, reason="bowtie2 not installed")
@pytest.mark.parametrize("gene,endogenous", [("renilla", False), ("aad4", True)])
def test_align_probes_bowtie2(gene, endogenous):
    """Verify bowtie2 produces valid SAM output."""
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

    if os.path.exists(fname_sam):
        os.remove(fname_sam)


# ===========================================================================
# 5. Masked intervals (original tests preserved)
# ===========================================================================

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
        assert task._is_in_masked_region("chrI", 150, 20) is True
        assert task._is_in_masked_region("chrI", 190, 20) is True
        assert task._is_in_masked_region("chrI", 300, 20) is False
        assert task._is_in_masked_region("chrIII", 150, 20) is False


# ###########################################################################
#  NEW UNIT TESTS -- all mocked, no external tools required
# ###########################################################################


class TestBowtie2CommandConstruction:
    """Verify the bowtie2 command differs for endogenous vs exogenous."""

    def _capture_bowtie2_args(self, is_endogenous, max_off_targets=0, min_length=25):
        """Run align_probes_bowtie2 with subprocess mocked; return the args list."""
        task = _make_task(
            is_endogenous=is_endogenous,
            max_off_targets=max_off_targets,
            aligner="bowtie2",
        )
        task.fname_genome = "/tmp/genome"

        # Create a small temp FASTA so Bio.SeqIO.parse works
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">probe1\nATGCATGCATGCATGCATGCA\n")
            task.fname_fasta = f.name

        task.fname_sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False).name

        captured = {}
        with patch("eFISHent.alignment.subprocess.check_call") as mock_cc, \
             patch("eFISHent.alignment.ProbeConfig") as mock_pc, \
             patch("eFISHent.alignment.os.remove"):
            mock_pc.return_value.min_length = min_length
            mock_cc.return_value = 0
            task.align_probes_bowtie2(threads=4)
            captured["args"] = mock_cc.call_args[0][0]

        os.unlink(task.fname_fasta)
        return captured["args"]

    def test_endogenous_uses_local_mode(self):
        args = self._capture_bowtie2_args(is_endogenous=True)
        assert "--local" in args
        assert "--end-to-end" not in args

    def test_exogenous_uses_end_to_end_mode(self):
        args = self._capture_bowtie2_args(is_endogenous=False)
        assert "--end-to-end" in args
        assert "--local" not in args

    def test_endogenous_score_min_log_based(self):
        args = self._capture_bowtie2_args(is_endogenous=True)
        idx = args.index("--score-min")
        assert args[idx + 1] == "G,1,4"

    def test_exogenous_score_min_linear(self):
        args = self._capture_bowtie2_args(is_endogenous=False)
        idx = args.index("--score-min")
        assert args[idx + 1] == "L,-0.6,-0.6"

    def test_k_value_minimum_100(self):
        """k should be at least 100 regardless of max_off_targets."""
        args = self._capture_bowtie2_args(is_endogenous=True, max_off_targets=0)
        idx = args.index("-k")
        assert int(args[idx + 1]) >= 100

    def test_k_value_scales_with_off_targets(self):
        """k = max(max_off_targets + 2, 100)."""
        args = self._capture_bowtie2_args(is_endogenous=True, max_off_targets=200)
        idx = args.index("-k")
        assert args[idx + 1] == "202"

    def test_fasta_input_flag(self):
        args = self._capture_bowtie2_args(is_endogenous=True)
        assert "-f" in args

    def test_threads_passed(self):
        args = self._capture_bowtie2_args(is_endogenous=True)
        idx = args.index("--threads")
        assert args[idx + 1] == "4"

    def test_adaptive_seed_length_short_probes(self):
        """Short probes (<=22) should use seed length 10."""
        args = self._capture_bowtie2_args(is_endogenous=True, min_length=20)
        idx = args.index("-L")
        assert args[idx + 1] == "10"

    def test_adaptive_seed_length_normal_probes(self):
        """Normal probes (>22) should use seed length 20."""
        args = self._capture_bowtie2_args(is_endogenous=True, min_length=25)
        idx = args.index("-L")
        assert args[idx + 1] == "20"


class TestBowtieCommandConstruction:
    """Verify bowtie (v1) command construction."""

    def _capture_bowtie_args(self, is_endogenous, max_off_targets=0):
        task = _make_task(
            is_endogenous=is_endogenous,
            max_off_targets=max_off_targets,
            aligner="bowtie",
        )
        task.fname_genome = "/tmp/genome"

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">probe1\nATGCATGCATGCATGCATGCA\n")
            task.fname_fasta = f.name

        task.fname_sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False).name

        with patch("eFISHent.alignment.subprocess.check_call") as mock_cc:
            mock_cc.return_value = 0
            task.align_probes(threads=4)
            args = mock_cc.call_args[0][0]

        os.unlink(task.fname_fasta)
        # Clean up the .fastq too
        fastq = task.fname_fasta.rstrip("a") + "q"
        if os.path.exists(fastq):
            os.unlink(fastq)
        return args

    def test_exogenous_has_m_flag(self):
        args = self._capture_bowtie_args(is_endogenous=False, max_off_targets=3)
        assert "-m" in args
        idx = args.index("-m")
        assert args[idx + 1] == "4"  # max_off_targets + 1

    def test_endogenous_no_m_flag(self):
        args = self._capture_bowtie_args(is_endogenous=True)
        assert "-m" not in args

    def test_all_flag_present(self):
        args = self._capture_bowtie_args(is_endogenous=True)
        assert "--all" in args

    def test_sam_flag_present(self):
        args = self._capture_bowtie_args(is_endogenous=True)
        assert "--sam" in args


class TestFilterUniqueProbesMocked:
    """Test filter_unique_probes with pysam.view mocked."""

    def _run_filter(self, sam_text, is_endogenous=True, aligner="bowtie2",
                    max_off_targets=0, no_alternative_loci=False):
        task = _make_task(
            is_endogenous=is_endogenous,
            aligner=aligner,
            max_off_targets=max_off_targets,
            no_alternative_loci=no_alternative_loci,
        )
        task.fname_sam = "/tmp/fake.sam"

        # Exogenous + bowtie2 tries to parse fname_fasta for the warning check
        tmpfasta = None
        if not is_endogenous and aligner == "bowtie2":
            tmpfasta = tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            )
            # Write enough probes so the warning logic can count them
            tmpfasta.write(">dummy\nATGC\n")
            tmpfasta.close()
            task.fname_fasta = tmpfasta.name

        with patch("eFISHent.alignment.pysam.view", return_value=sam_text):
            result = task.filter_unique_probes()

        if tmpfasta is not None:
            os.unlink(tmpfasta.name)

        return result

    def test_endogenous_bowtie2_excludes_unmapped(self):
        """Endogenous + bowtie2: pysam.view called with --exclude-flags 4."""
        task = _make_task(is_endogenous=True, aligner="bowtie2")
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value="") as mock_view:
            task.filter_unique_probes()
            mock_view.assert_called_once_with("/tmp/fake.sam", "--exclude-flags", "4")

    def test_exogenous_bowtie2_requires_unmapped(self):
        """Exogenous + bowtie2: pysam.view called with --require-flags 4."""
        task = _make_task(is_endogenous=False, aligner="bowtie2")
        task.fname_sam = "/tmp/fake.sam"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">probe1\nATGC\n")
            task.fname_fasta = f.name

        with patch("eFISHent.alignment.pysam.view", return_value="") as mock_view:
            task.filter_unique_probes()
            mock_view.assert_called_once_with("/tmp/fake.sam", "--require-flags", "4")

        os.unlink(task.fname_fasta)

    def test_endogenous_bowtie_uses_mapq60(self):
        """Endogenous + bowtie: pysam.view called with --min-MQ 60."""
        task = _make_task(is_endogenous=True, aligner="bowtie")
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value="") as mock_view:
            task.filter_unique_probes()
            mock_view.assert_called_once_with("/tmp/fake.sam", "--min-MQ", "60")

    def test_exogenous_bowtie_requires_unmapped(self):
        """Exogenous + bowtie: pysam.view called with --require-flags 4."""
        task = _make_task(is_endogenous=False, aligner="bowtie")
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value="") as mock_view:
            task.filter_unique_probes()
            mock_view.assert_called_once_with("/tmp/fake.sam", "--require-flags", "4")

    def test_no_alignments_returns_empty(self):
        df, df_pre = self._run_filter("")
        assert len(df) == 0
        assert len(df_pre) == 0

    def test_single_unique_probe_kept(self):
        """A probe with exactly one hit should be kept (endogenous, max_off_targets=0)."""
        sam = _make_sam_line("probe1", "0", "chrI", "1000", "60")
        df, _ = self._run_filter(sam, is_endogenous=True, max_off_targets=0)
        assert len(df) == 1
        assert df.iloc[0]["qname"] == "probe1"

    def test_off_target_filtering_endogenous(self):
        """Probe with too many hits should be removed in endogenous mode."""
        # probe1 has 3 hits, max_off_targets=1 means keep <= 2 hits
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
            _make_sam_line("probe1", "0", "chrIII", "3000", "60"),
        ])
        df, _ = self._run_filter(lines, is_endogenous=True, max_off_targets=1)
        # 3 hits > max_off_targets+1 (2), so probe1 should be filtered out
        assert len(df) == 0

    def test_off_target_filtering_keeps_within_threshold(self):
        """Probe with hits within threshold should be kept."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
        ])
        df, _ = self._run_filter(lines, is_endogenous=True, max_off_targets=1)
        # 2 hits <= max_off_targets+1 (2), so probe1 should be kept
        assert len(df) == 2
        assert df.iloc[0]["qname"] == "probe1"

    def test_mixed_probes_some_filtered(self):
        """Some probes pass, some fail off-target threshold."""
        lines = "\n".join([
            # probe1: 1 hit -> kept
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            # probe2: 3 hits -> removed (max_off_targets=0, threshold=1)
            _make_sam_line("probe2", "0", "chrI", "2000", "60"),
            _make_sam_line("probe2", "0", "chrII", "3000", "60"),
            _make_sam_line("probe2", "0", "chrIII", "4000", "60"),
        ])
        df, df_pre = self._run_filter(lines, is_endogenous=True, max_off_targets=0)
        assert set(df["qname"].unique()) == {"probe1"}
        # df_pre_count should still have all probes (snapshot before count filter)
        assert set(df_pre["qname"].unique()) == {"probe1", "probe2"}

    def test_alternative_loci_filtered(self):
        """Alignments to _alt contigs should be removed when no_alternative_loci=True."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrI_alt", "2000", "60"),
        ])
        df, _ = self._run_filter(
            lines, is_endogenous=True, max_off_targets=0, no_alternative_loci=True
        )
        # After removing _alt hit, probe1 has 1 hit <= threshold 1
        assert len(df) == 1
        assert "chrI_alt" not in df["rname"].values

    def test_alternative_loci_not_filtered_when_disabled(self):
        """_alt contigs kept when no_alternative_loci=False."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrI_alt", "2000", "60"),
        ])
        df, _ = self._run_filter(
            lines, is_endogenous=True, max_off_targets=1, no_alternative_loci=False
        )
        assert len(df) == 2

    def test_exogenous_mode_no_count_filtering(self):
        """Exogenous mode does not apply endogenous count-based filtering."""
        lines = "\n".join([
            _make_sam_line("probe1", "4", "*", "0", "0"),
            _make_sam_line("probe2", "4", "*", "0", "0"),
        ])
        df, _ = self._run_filter(lines, is_endogenous=False, max_off_targets=0)
        assert len(df) == 2

    def test_pre_count_snapshot_includes_all(self):
        """df_pre_count should include everything before count-based filtering."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe2", "0", "chrI", "2000", "60"),
            _make_sam_line("probe2", "0", "chrII", "3000", "60"),
            _make_sam_line("probe2", "0", "chrIII", "4000", "60"),
        ])
        df, df_pre = self._run_filter(lines, is_endogenous=True, max_off_targets=0)
        assert set(df_pre["qname"].unique()) == {"probe1", "probe2"}
        assert set(df["qname"].unique()) == {"probe1"}

    def test_bowtie2_partial_alignments_filtered_by_cigar_match_fraction(self):
        """Short local alignments should not count as real off-targets."""
        lines = "\n".join([
            "probe1\t0\tchrI\t1000\t60\t21M\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*",
            "probe1\t0\tchrII\t2000\t60\t10M11S\t*\t0\t0\tATGCATGCATGCATGCATGCA\t*",
        ])
        df, df_pre = self._run_filter(lines, is_endogenous=True, max_off_targets=0)
        assert set(df_pre["cigar"]) == {"21M"}
        assert len(df) == 1
        assert df.iloc[0]["rname"] == "chrI"

    def test_transcriptome_presence_skips_genome_count_filter_at_default(self):
        """With transcriptome BLAST enabled, default genome hit count should not reject probes."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
            _make_sam_line("probe1", "0", "chrIII", "3000", "60"),
        ])
        task = _make_task(is_endogenous=True, aligner="bowtie2", max_off_targets=0)
        task.fname_sam = "/tmp/fake.sam"
        task.has_transcriptome = True
        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, df_pre = task.filter_unique_probes()
        assert len(df_pre) == 3
        assert len(df) == 3


class TestAlignmentHelpers:
    """Unit tests for helper filters added around alignment parsing."""

    def test_cigar_match_fraction_counts_soft_clips_and_insertions_in_denominator(self):
        assert AlignProbeCandidates._cigar_match_fraction("16M2I3S") == pytest.approx(16 / 21)

    def test_filter_pseudogene_off_targets_removes_overlapping_hits(self, tmp_path):
        task = _make_task(is_endogenous=True, aligner="bowtie2")
        task.has_pseudogene_filter = True

        annotation_path = tmp_path / "annotation.gtf.parq"
        pd.DataFrame(
            {
                "seqname": ["chr1", "chr1", "chr2"],
                "start": [90, 500, 700],
                "end": [150, 550, 750],
                "feature": ["gene", "gene", "gene"],
                "gene_biotype": ["processed_pseudogene", "protein_coding", "pseudogene"],
            }
        ).to_parquet(annotation_path)

        task.input = lambda: {"annotation": MagicMock(path=str(annotation_path))}
        df = pd.DataFrame(
            {
                "qname": ["probe1", "probe1", "probe2"],
                "rname": ["chr1", "chr1", "chr2"],
                "pos": ["100", "520", "720"],
                "cigar": ["21M", "21M", "21M"],
            }
        )

        filtered = task.filter_pseudogene_off_targets(df)
        assert len(filtered) == 1
        assert filtered.iloc[0]["pos"] == "520"


class TestCountTableParsing:
    """Test count_table property with various file formats."""

    def test_csv_extension(self):
        """CSV files should be parsed with comma separator."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            # index_col=0 consumes the first column as index, so after
            # reset_index + iloc[:, :2] the first two are gene_id, count
            f.write("gene_id,count\nENSG00000001.1,42.5\nENSG00000002.3,10.0\n")
            task.fname_count = f.name

        df = task.count_table
        os.unlink(task.fname_count)

        assert list(df.columns) == ["gene_id", "count"]
        # Version numbers should be stripped
        assert "ENSG00000001" in df["gene_id"].values
        assert "ENSG00000002" in df["gene_id"].values

    def test_tsv_extension(self):
        """TSV files should be parsed with tab separator."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("gene_id\tcount\nENSG00000001\t42.5\n")
            task.fname_count = f.name

        df = task.count_table
        os.unlink(task.fname_count)

        assert len(df) == 1
        assert df.iloc[0]["count"] == 42.5

    def test_invalid_extension_raises(self):
        """Non-TSV/CSV/TXT files should raise ValueError."""
        task = _make_task()
        task.fname_count = "/tmp/data.json"
        with pytest.raises(ValueError, match="TSV, CSV, or TXT"):
            _ = task.count_table

    def test_version_numbers_stripped(self):
        """Gene ID version suffixes (e.g., .3) should be removed."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("gene_id,count\nENSG00000139618.15,100\n")
            task.fname_count = f.name

        df = task.count_table
        os.unlink(task.fname_count)

        assert df.iloc[0]["gene_id"] == "ENSG00000139618"

    def test_single_column_raises(self):
        """Count table with only one column should raise ValueError."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("only_col\n42\n")
            task.fname_count = f.name

        with pytest.raises(ValueError, match="at least 2 columns"):
            _ = task.count_table
        os.unlink(task.fname_count)

    def test_non_numeric_count_raises(self):
        """Non-numeric count column should raise ValueError."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("gene_id,count\nGENE1,notanumber\n")
            task.fname_count = f.name

        with pytest.raises(ValueError, match="convert count table"):
            _ = task.count_table
        os.unlink(task.fname_count)

    def test_na_rows_dropped(self):
        """Rows with NaN values should be dropped."""
        task = _make_task()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("gene_id,count\nGENE1,10\n,\nGENE3,30\n")
            task.fname_count = f.name

        df = task.count_table
        os.unlink(task.fname_count)

        assert len(df) == 2


class TestFilterMaskedOffTargets:
    """Test filter_masked_off_targets with dustmasker mocked."""

    def test_removes_masked_hits(self):
        """Hits in masked regions should be removed."""
        task = _make_task(is_endogenous=True, mask_repeats=True)

        df = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "150", "seq": "ATGCATGCATGCATGCATGCA"},
            {"qname": "probe1", "rname": "chrI", "pos": "800", "seq": "ATGCATGCATGCATGCATGCA"},
        ])

        task._masked_intervals = {"chrI": ([100], [200])}

        with patch.object(task, "_run_dustmasker", return_value="/tmp/mask.intervals"), \
             patch.object(
                 AlignProbeCandidates, "_load_masked_intervals",
                 return_value={"chrI": ([100], [200])}
             ), \
             patch("eFISHent.alignment.shutil.which", return_value="/usr/bin/dustmasker"):
            result = task.filter_masked_off_targets(df)

        assert len(result) == 1
        assert result.iloc[0]["pos"] == "800"

    def test_skips_when_dustmasker_missing(self):
        """Should return unmodified df when dustmasker is not installed."""
        task = _make_task(mask_repeats=True)
        df = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "150", "seq": "ATGC"},
        ])

        with patch("eFISHent.alignment.shutil.which", return_value=None):
            result = task.filter_masked_off_targets(df)

        assert len(result) == 1


class TestFilterIntergenicOffTargets:
    """Test filter_intergenic_off_targets with GTF mocked."""

    def test_removes_intergenic_hits(self):
        """Hits not overlapping any gene should be removed."""
        task = _make_task(is_endogenous=True, has_intergenic_filter=True)

        gtf_df = pd.DataFrame({
            "gene_id": ["GENE1"],
            "seqname": ["chrI"],
            "start": [1000],
            "end": [2000],
            "frame": ["."],
            "exon_number": [1],
            "feature": ["gene"],
        })

        mock_input = MagicMock()
        mock_input.__getitem__ = MagicMock(return_value=MagicMock(path="/tmp/ann.parquet"))
        task.input = MagicMock(return_value=mock_input)

        df_sam = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "1500"},   # inside gene
            {"qname": "probe2", "rname": "chrI", "pos": "5000"},   # intergenic
        ])

        with patch("eFISHent.alignment.pd.read_parquet", return_value=gtf_df):
            result = task.filter_intergenic_off_targets(df_sam)

        assert len(result) == 1
        assert result.iloc[0]["qname"] == "probe1"

    def test_no_genes_removes_all(self):
        """If no genes annotated, all hits should be removed."""
        task = _make_task(has_intergenic_filter=True)

        gtf_df = pd.DataFrame(
            columns=["gene_id", "seqname", "start", "end", "frame", "exon_number", "feature"]
        )

        mock_input = MagicMock()
        mock_input.__getitem__ = MagicMock(return_value=MagicMock(path="/tmp/ann.parquet"))
        task.input = MagicMock(return_value=mock_input)

        df_sam = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "1500"},
        ])

        with patch("eFISHent.alignment.pd.read_parquet", return_value=gtf_df):
            result = task.filter_intergenic_off_targets(df_sam)

        assert len(result) == 0


class TestFilterRrnaOffTargets:
    """Test filter_rrna_off_targets with annotation mocked."""

    def _build_gtf(self, biotype_col="gene_biotype"):
        return pd.DataFrame({
            "gene_id": ["RRNA1", "GENE2"],
            biotype_col: ["rRNA", "protein_coding"],
            "feature": ["gene", "gene"],
            "seqname": ["chrI", "chrII"],
            "start": [100, 5000],
            "end": [500, 6000],
            "frame": [".", "."],
            "exon_number": [1, 1],
        })

    def test_removes_probes_hitting_rrna(self):
        """Probes with any rRNA hit should be entirely removed."""
        task = _make_task(has_rrna_filter=True, filter_rrna=True)
        gtf_df = self._build_gtf()

        mock_input = MagicMock()
        mock_input.__getitem__ = MagicMock(return_value=MagicMock(path="/tmp/ann.parquet"))
        task.input = MagicMock(return_value=mock_input)

        df_sam = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "200"},   # hits rRNA
            {"qname": "probe1", "rname": "chrII", "pos": "5500"},  # also maps elsewhere
            {"qname": "probe2", "rname": "chrII", "pos": "5500"},  # safe
        ])

        with patch("eFISHent.alignment.pd.read_parquet", return_value=gtf_df):
            result = task.filter_rrna_off_targets(df_sam)

        # probe1 entirely removed (both rows), probe2 kept
        assert set(result["qname"].unique()) == {"probe2"}

    def test_no_biotype_column_skips(self):
        """If GTF has no biotype column, skip filter gracefully."""
        task = _make_task(has_rrna_filter=True)
        gtf_df = pd.DataFrame({
            "gene_id": ["GENE1"],
            "feature": ["gene"],
            "seqname": ["chrI"],
            "start": [100],
            "end": [500],
        })

        mock_input = MagicMock()
        mock_input.__getitem__ = MagicMock(return_value=MagicMock(path="/tmp/ann.parquet"))
        task.input = MagicMock(return_value=mock_input)

        df_sam = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "200"},
        ])

        with patch("eFISHent.alignment.pd.read_parquet", return_value=gtf_df):
            result = task.filter_rrna_off_targets(df_sam)

        assert len(result) == 1

    def test_no_rrna_genes_keeps_all(self):
        """If no rRNA genes exist, all probes should be kept."""
        task = _make_task(has_rrna_filter=True)
        gtf_df = pd.DataFrame({
            "gene_id": ["GENE1"],
            "gene_biotype": ["protein_coding"],
            "feature": ["gene"],
            "seqname": ["chrI"],
            "start": [100],
            "end": [500],
            "frame": ["."],
            "exon_number": [1],
        })

        mock_input = MagicMock()
        mock_input.__getitem__ = MagicMock(return_value=MagicMock(path="/tmp/ann.parquet"))
        task.input = MagicMock(return_value=mock_input)

        df_sam = pd.DataFrame([
            {"qname": "probe1", "rname": "chrI", "pos": "200"},
        ])

        with patch("eFISHent.alignment.pd.read_parquet", return_value=gtf_df):
            result = task.filter_rrna_off_targets(df_sam)

        assert len(result) == 1


class TestMaskedIntervalsExtended:
    """Additional edge-case tests for masked interval logic."""

    def test_empty_interval_file(self):
        """Empty interval file should yield empty dict."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".intervals", delete=False) as f:
            f.write("")
            path = f.name

        result = AlignProbeCandidates._load_masked_intervals(path)
        os.unlink(path)
        assert result == {}

    def test_header_only_no_intervals(self):
        """Chromosome header with no intervals should not appear."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".intervals", delete=False) as f:
            f.write(">chrI\n>chrII\n100 - 200\n")
            path = f.name

        result = AlignProbeCandidates._load_masked_intervals(path)
        os.unlink(path)

        # chrI had no intervals before chrII header
        assert "chrI" not in result
        assert "chrII" in result

    def test_boundary_overlap_start(self):
        """Read starting just before masked region end should overlap."""
        task = _make_task()
        task._masked_intervals = {"chrI": ([100], [200])}
        # Read at pos=195, length=20 -> [195, 215) overlaps [100, 200]
        assert task._is_in_masked_region("chrI", 195, 20) is True

    def test_boundary_no_overlap(self):
        """Read starting after masked region should not overlap."""
        task = _make_task()
        task._masked_intervals = {"chrI": ([100], [200])}
        # Read at pos=201, length=20 -> [201, 221) does not overlap [100, 200]
        assert task._is_in_masked_region("chrI", 201, 20) is False


class TestEdgeCases:
    """Edge cases for the alignment pipeline."""

    def test_all_probes_filtered_endogenous(self):
        """When all probes exceed off-target threshold, result should be empty."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
            _make_sam_line("probe2", "0", "chrI", "3000", "60"),
            _make_sam_line("probe2", "0", "chrIII", "4000", "60"),
        ])
        task = _make_task(is_endogenous=True, max_off_targets=0)
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, _ = task.filter_unique_probes()

        assert len(df) == 0

    def test_single_hit_per_probe_all_pass(self):
        """Each probe with exactly one hit should pass max_off_targets=0."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe2", "0", "chrII", "2000", "60"),
            _make_sam_line("probe3", "0", "chrIII", "3000", "60"),
        ])
        task = _make_task(is_endogenous=True, max_off_targets=0)
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, _ = task.filter_unique_probes()

        assert len(df) == 3
        assert set(df["qname"].unique()) == {"probe1", "probe2", "probe3"}

    def test_exogenous_all_unmapped_warning(self):
        """Exogenous mode should warn when all probes pass as unmapped."""
        lines = "\n".join([
            _make_sam_line("probe1", "4", "*", "0", "0"),
            _make_sam_line("probe2", "4", "*", "0", "0"),
        ])
        task = _make_task(is_endogenous=False, aligner="bowtie2")
        task.fname_sam = "/tmp/fake.sam"

        # Create a temp FASTA with 2 probes
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(">probe1\nATGC\n>probe2\nGCTA\n")
            task.fname_fasta = f.name

        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, _ = task.filter_unique_probes()

        os.unlink(task.fname_fasta)
        # All probes should still pass
        assert len(df) == 2

    def test_max_off_targets_boundary(self):
        """Probe at exactly the threshold should be kept."""
        # max_off_targets=2 means keep probes with <= 3 hits (2 off-target + 1 on-target)
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
            _make_sam_line("probe1", "0", "chrIII", "3000", "60"),
        ])
        task = _make_task(is_endogenous=True, max_off_targets=2)
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, _ = task.filter_unique_probes()

        # 3 hits <= max_off_targets+1 (3) -> kept
        assert len(df) == 3

    def test_max_off_targets_one_over_boundary(self):
        """Probe with one hit over threshold should be removed."""
        lines = "\n".join([
            _make_sam_line("probe1", "0", "chrI", "1000", "60"),
            _make_sam_line("probe1", "0", "chrII", "2000", "60"),
            _make_sam_line("probe1", "0", "chrIII", "3000", "60"),
            _make_sam_line("probe1", "0", "chrIV", "4000", "60"),
        ])
        task = _make_task(is_endogenous=True, max_off_targets=2)
        task.fname_sam = "/tmp/fake.sam"

        with patch("eFISHent.alignment.pysam.view", return_value=lines):
            df, _ = task.filter_unique_probes()

        # 4 hits > max_off_targets+1 (3) -> removed
        assert len(df) == 0
