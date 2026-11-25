"""Tests for the session cleanup in tests/conftest.py."""

from tests.conftest import remove_stale_test_files


def test_cleanup_removes_only_generated_indices(tmp_path):
    kept = ["sacCer3.fa", "sacCer3.gtf", "sacCer3.gtf.parq", "renilla_basic.fasta"]
    stale = ["sacCer3.1.ebwt", "sacCer3.rev.1.bt2", "sacCer3_15.jf"]
    for name in kept + stale:
        (tmp_path / name).write_text("x")

    remove_stale_test_files(str(tmp_path))

    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(kept)
