"""Shared pytest hooks for the eFISHent test suite."""

import glob
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Indices that test_integration.py builds next to tests/data/sacCer3.fa
# (bowtie, bowtie2, jellyfish), deleted so a rerun builds them fresh.
# The tracked inputs (sacCer3.fa, sacCer3.gtf, sacCer3.gtf.parq, ...) never match.
STALE_PATTERNS = ("sacCer3*.ebwt", "sacCer3*.bt2", "sacCer3*.jf")


def remove_stale_test_files(data_dir: str = DATA_DIR) -> None:
    """Delete indices left in data_dir by an earlier run, so tests rebuild them."""
    for pattern in STALE_PATTERNS:
        for path in glob.glob(os.path.join(data_dir, pattern)):
            if os.path.isfile(path):
                os.remove(path)


# Runs once after collection, before any test. Only when integration tests are
# selected: they rebuild the indices, while unit tests (kmers, alignment,
# cleanup) read the indices an earlier integration run left behind.
def pytest_collection_finish(session):
    if any(item.path.name == "test_integration.py" for item in session.items):
        remove_stale_test_files()
