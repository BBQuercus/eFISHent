import os
import shutil
import sys
from unittest.mock import patch

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.secondary_structure import get_free_energy

# Check if Fold is available (bundled binary or on PATH)
_fold_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_fold_bin = os.path.join(_fold_dir, "eFISHent", "Fold_linux" if sys.platform.startswith("linux") else "Fold_osx")
FOLD_AVAILABLE = os.access(_fold_bin, os.X_OK) or shutil.which("Fold") is not None


@pytest.mark.skipif(not FOLD_AVAILABLE, reason="Fold binary not available")
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


def test_get_free_energy_fold_failure_returns_zero():
    """When Fold binary exits with non-zero status, return 0.0 instead of crashing."""
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="test")

    with patch("eFISHent.secondary_structure.subprocess.run") as mock_run:
        import subprocess
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd=["Fold", "-", "-"], stderr="some error"
        )
        energy = get_free_energy(sequence)
    assert energy == 0.0


def test_get_free_energy_fold_failure_logs_warning():
    """When Fold binary fails, a warning should be logged."""
    sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="test_seq")

    with patch("eFISHent.secondary_structure.subprocess.run") as mock_run:
        import subprocess
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd=["Fold", "-", "-"], stderr="segfault"
        )
        import logging
        with patch.object(logging.getLogger("custom-logger"), "warning") as mock_warn:
            get_free_energy(sequence)
            mock_warn.assert_called_once()
            assert "test_seq" in mock_warn.call_args[0][0]


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
