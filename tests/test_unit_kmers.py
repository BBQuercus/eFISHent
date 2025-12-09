import shutil

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.kmers import get_max_kmer_count
from eFISHent.kmers import get_max_kmer_counts_batch

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
