from unittest.mock import patch

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


def test_get_free_energy_unsupported_platform():
    """Test that unsupported platforms raise NotImplementedError.

    This test exposes the bug where Windows (and other platforms) would
    cause a NameError because fold_path is never defined.
    The fix should raise a clear NotImplementedError instead.
    """
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="test")

    with patch("eFISHent.secondary_structure.sys.platform", "win32"):
        with pytest.raises(NotImplementedError, match="not supported"):
            get_free_energy(sequence)
