import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.secondary_structure import get_free_energy


@pytest.mark.parametrize(
    "seq,deltag",
    [
        ("AACTTGTCTTAGCTTTGCAGTCGAGTT", -4.8),
        ("TAGCTTTGC", 0.0),
        ("ACGTGCCACGATTCAACGTGGCACAG", -15.1),
    ],
)
def test_get_free_energy(seq, deltag):
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id="sequence")
    assert get_free_energy(sequence) == deltag
