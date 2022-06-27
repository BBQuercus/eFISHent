"""Unique and verified Length/TM/GC/formamid/na/off-targets/kmers/deltaG."""

import glob
import os
import subprocess

import pytest


def verify_csv(fname: str, *params):
    pass


def test_passed_exogenous():
    """Exogenous sequence (renilla) as sequence file, unique output dir, all files saved."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--save-intermediates",
        "True",
        "--sequence-file",
        "./tests/renilla.fa",
        "--output-dir",
        "./tests/",
        "--is-endogenous",
        "False",
        "--min-length",
        "20",
        "--max-length",
        "24",
        "--min-tm",
        "45",
        "--max-tm",
        "55",
    ]
    process = subprocess.run(args)
    assert process.returncode == 0
    files = glob.glob("./renilla_*.csv")
    assert files
    [os.remove(f) for f in files]


def test_downloaded_endogenous_optimal():
    """Gene/organism downloaded downloaded, minus strand, optimal model w/ timeout."""
    pass


def test_counttable_full_output():
    """Ensembl downloaded sequence with count table and max off target FPKM."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--ensembl-id",
        "YKL185W",
    ]
    process = subprocess.run(args)
    assert process.returncode == 0
    files = glob.glob("./YKL185W_*")
    assert len(files) == 2
    [os.remove(f) for f in files]


def test_index_only():
    """Build index step only."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--build-indices",
        "True",
    ]
    process = subprocess.run(args)
    assert process.returncode == 0


def test_help_only():
    """Help only with and without being explicit."""
    process_1 = subprocess.check_output("efishent").decode()
    process_2 = subprocess.check_output("eFISHent").decode()
    process_3 = subprocess.check_output(["efishent", "--help"]).decode()
    assert process_1 == process_2 == process_3


def test_error():
    """No reference genome passed, bad values (too short length, non-bools)."""
    # Wrong argument names
    process = subprocess.run(["eFISHent", "--not-an-argument", "ENSB29329382938"])
    assert process.returncode == 2

    # No reference genome passed
    process = subprocess.run(["eFISHent", "--ensembl-id", "ENSG00000026025"])
    assert process.returncode == 2
    files = glob.glob("./ENSG00000026025_*.csv")
    assert len(files) == 0

    # Non-bools
    process = subprocess.run(
        [
            "efishent",
            "--reference-genome",
            "./sacCer3.fa",
            "--save-intermediates",
            "Yup",
            "--sequence-file",
            "./renilla.fa",
        ]
    )
    assert process.returncode == 2

    # Too short/long length?
    # p = subprocess.Popen(args, stderr=subprocess.PIPE)
    # stdout, stderr = p.communicate()
