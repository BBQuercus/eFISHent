"""Unique and verified Length/TM/GC/formamid/na/off-targets/kmers/deltaG."""

import subprocess

import pytest


def verify_csv(fname: str, *params):
    pass


def test_passed_exogenous():
    """Exogenous sequence (renilla) as sequence file, unique output dir, all files saved."""
    args = [
        "eFISHent",
        "--reference-genome",
        "./sacCer3.fa",
        "--save-intermediates",
        "True",
        "--sequence-file",
        "./renilla.fa",
        "--output-dir",
        "./?",
        "--min-length",
        "20",
        "--max-length",
        "24",
        "--min-tm",
        "45",
        "--max-tm",
        "55",
    ]
    # subprocess.check_call(args)


def test_downloaded_endogenous_optimal():
    """Gene/organism downloaded downloaded, minus strand, optimal model w/ timeout."""
    pass


def test_counttable_full_output():
    """Ensembl downloaded sequence with count table and max off target FPKM."""
    pass


def test_index_only():
    """Build index step only."""
    pass


def test_help_only():
    """Help only with and without being explicit."""
    pass


def test_error():
    """No reference genome passed, bad values (too short length, non-bools)."""
    pass
