import shutil

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.kmers import get_max_kmer_count

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
