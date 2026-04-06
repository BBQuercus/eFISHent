import shutil
from unittest.mock import MagicMock, patch

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.kmers import KMerFiltering, get_max_kmer_count, get_max_kmer_counts_batch

JELLYFISH_AVAILABLE = shutil.which("jellyfish") is not None


@pytest.mark.skipif(not JELLYFISH_AVAILABLE, reason="jellyfish not installed")
@pytest.mark.parametrize(
    "seq,count",
    [
        ("ATATATATATATATATATATATAT", 255),
        ("GGGGGGGGGGGGGGGGGGGGGGGGG", 2),
        ("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 2),
        ("ATCACTAGTACGTAGATGTGCG", 0),
    ],
)
def test_max_kmer_count(seq, count):
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id="seq")
    assert get_max_kmer_count(sequence, "./tests/data/sacCer3_15.jf") == count


@pytest.mark.parametrize("seq", ["ATAT", "GCATGCATGC"])
def test_max_kmer_count_error(seq):
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id="seq")
    with pytest.raises(ValueError):
        get_max_kmer_count(sequence, "./tests/data/sacCer3_15.jf")


@pytest.mark.skipif(not JELLYFISH_AVAILABLE, reason="jellyfish not installed")
def test_max_kmer_counts_batch():
    """Test batch kmer counting returns same results as individual calls.

    This test ensures the batched implementation produces identical results
    to individual calls while being more efficient (single subprocess call).
    """
    sequences = [
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATATATATATATATATATATATAT"), id="seq1"),
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("GGGGGGGGGGGGGGGGGGGGGGGGG"), id="seq2"),
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCACTAGTACGTAGATGTGCG"), id="seq3"),
    ]
    expected_counts = [255, 2, 0]

    # Batch function should return same results as individual calls
    batch_counts = get_max_kmer_counts_batch(sequences, "./tests/data/sacCer3_15.jf")

    assert batch_counts == expected_counts


@pytest.mark.skipif(not JELLYFISH_AVAILABLE, reason="jellyfish not installed")
def test_max_kmer_counts_batch_empty():
    """Test batch kmer counting with empty list returns empty list."""
    batch_counts = get_max_kmer_counts_batch([], "./tests/data/sacCer3_15.jf")
    assert batch_counts == []


# ── Mocked get_max_kmer_count tests ─────────────────────────────────────


def _make_seq(seq_str, seq_id="seq"):
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq_str), id=seq_id)


class TestGetMaxKmerCountMocked:
    """Test get_max_kmer_count with mocked subprocess (no jellyfish needed)."""

    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.subprocess.check_output")
    def test_correct_max_extraction(self, mock_subprocess, mock_config):
        mock_config.return_value.kmer_length = 5
        # Simulate jellyfish output: 3 kmers with counts 3, 7, 2
        mock_subprocess.return_value = b"ATGCA 3\nTGCAT 7\nGCATG 2\n"

        seq = _make_seq("ATGCATGCA")  # length 9, kmer_length 5 -> 5 kmers
        result = get_max_kmer_count(seq, "/fake/path.jf")
        assert result == 7

    @patch("eFISHent.kmers.ProbeConfig")
    def test_sequence_shorter_than_kmer(self, mock_config):
        mock_config.return_value.kmer_length = 15
        seq = _make_seq("ATGCATGC")  # length 8 < 15
        with pytest.raises(ValueError):
            get_max_kmer_count(seq, "/fake/path.jf")

    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.subprocess.check_output")
    def test_all_zeros(self, mock_subprocess, mock_config):
        mock_config.return_value.kmer_length = 4
        mock_subprocess.return_value = b"ATGC 0\nTGCA 0\nGCAT 0\n"

        seq = _make_seq("ATGCATG")  # length 7, kmer_length 4 -> 4 kmers
        result = get_max_kmer_count(seq, "/fake/path.jf")
        assert result == 0


# ── Mocked get_max_kmer_counts_batch tests ───────────────────────────────


class TestGetMaxKmerCountsBatchMocked:
    """Test get_max_kmer_counts_batch with mocked subprocess."""

    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.subprocess.check_output")
    def test_batch_matches_individual(self, mock_subprocess, mock_config):
        mock_config.return_value.kmer_length = 4
        # Two sequences: "ATGCATG" (4 kmers) and "GGCCGG" (3 kmers)
        # jellyfish output for all 7 kmers combined
        mock_subprocess.return_value = (
            b"ATGC 2\nTGCA 5\nGCAT 1\nCATG 3\n"  # seq1: max=5
            b"GGCC 10\nGCCG 4\nCCGG 8\n"  # seq2: max=10
        )

        sequences = [_make_seq("ATGCATG", "s1"), _make_seq("GGCCGG", "s2")]
        result = get_max_kmer_counts_batch(sequences, "/fake/path.jf")
        assert result == [5, 10]

    def test_empty_list(self):
        result = get_max_kmer_counts_batch([], "/fake/path.jf")
        assert result == []

    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.subprocess.check_output")
    def test_single_sequence_batch(self, mock_subprocess, mock_config):
        mock_config.return_value.kmer_length = 4
        mock_subprocess.return_value = b"ATGC 3\nTGCA 7\nGCAT 1\nCATG 2\n"

        sequences = [_make_seq("ATGCATG", "s1")]
        result = get_max_kmer_counts_batch(sequences, "/fake/path.jf")
        assert result == [7]


# ── KMerFiltering task logic (mocked) ────────────────────────────────────


class TestKMerFilteringLogic:
    """Test KMerFiltering run logic with mocked dependencies."""

    def _make_probes(self, n=5):
        """Create n dummy probe SeqRecords."""
        return [
            _make_seq("ATGCATGCATGCATGC", f"probe-{i}") for i in range(n)
        ]

    @patch("eFISHent.kmers.Bio.SeqIO")
    @patch("eFISHent.kmers.get_max_kmer_counts_batch")
    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.GeneralConfig")
    @patch("eFISHent.kmers.util")
    def test_exogenous_skips_filtering(
        self, mock_util, mock_general, mock_probe_config, mock_batch, mock_seqio
    ):
        """Exogenous genes should skip kmer filtering entirely."""
        probes = self._make_probes(5)
        mock_seqio.parse.return_value = probes

        task = KMerFiltering()
        task.input = MagicMock(return_value={
            "probes": {"fasta": MagicMock(path="/fake/probes.fasta")},
            "jellyfish": MagicMock(path="/fake/index.jf"),
        })
        task.output = MagicMock(return_value=MagicMock(path="/fake/output.fasta"))

        with patch("eFISHent.config.SequenceConfig") as mock_seq_config:
            mock_seq_config.return_value.is_endogenous = False
            task.run()

        # Batch kmer counting should NOT be called for exogenous
        mock_batch.assert_not_called()
        # All probes should be written out
        written_probes = mock_seqio.write.call_args[0][0]
        assert len(written_probes) == 5

    @patch("eFISHent.kmers.Bio.SeqIO")
    @patch("eFISHent.kmers.get_max_kmer_counts_batch")
    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.GeneralConfig")
    @patch("eFISHent.kmers.util")
    def test_endogenous_filters_by_threshold(
        self, mock_util, mock_general, mock_probe_config, mock_batch, mock_seqio
    ):
        """Endogenous genes should filter probes exceeding max_kmers."""
        probes = self._make_probes(4)
        mock_seqio.parse.return_value = probes
        mock_probe_config.return_value.max_kmers = 5

        # Counts: probe-0=3 (pass), probe-1=10 (fail), probe-2=5 (pass, at boundary), probe-3=6 (fail)
        mock_batch.return_value = [3, 10, 5, 6]

        task = KMerFiltering()
        task.input = MagicMock(return_value={
            "probes": {"fasta": MagicMock(path="/fake/probes.fasta")},
            "jellyfish": MagicMock(path="/fake/index.jf"),
        })
        task.output = MagicMock(return_value=MagicMock(path="/fake/output.fasta"))

        with patch("eFISHent.config.SequenceConfig") as mock_seq_config:
            mock_seq_config.return_value.is_endogenous = True
            task.run()

        written_probes = mock_seqio.write.call_args[0][0]
        assert len(written_probes) == 2
        assert written_probes[0].id == "probe-0"
        assert written_probes[1].id == "probe-2"

    @patch("eFISHent.kmers.Bio.SeqIO")
    @patch("eFISHent.kmers.get_max_kmer_counts_batch")
    @patch("eFISHent.kmers.ProbeConfig")
    @patch("eFISHent.kmers.GeneralConfig")
    @patch("eFISHent.kmers.util")
    def test_boundary_exactly_at_max_kmers(
        self, mock_util, mock_general, mock_probe_config, mock_batch, mock_seqio
    ):
        """Probes with count exactly equal to max_kmers should pass."""
        probes = self._make_probes(3)
        mock_seqio.parse.return_value = probes
        mock_probe_config.return_value.max_kmers = 4

        mock_batch.return_value = [4, 4, 4]

        task = KMerFiltering()
        task.input = MagicMock(return_value={
            "probes": {"fasta": MagicMock(path="/fake/probes.fasta")},
            "jellyfish": MagicMock(path="/fake/index.jf"),
        })
        task.output = MagicMock(return_value=MagicMock(path="/fake/output.fasta"))

        with patch("eFISHent.config.SequenceConfig") as mock_seq_config:
            mock_seq_config.return_value.is_endogenous = True
            task.run()

        written_probes = mock_seqio.write.call_args[0][0]
        assert len(written_probes) == 3
