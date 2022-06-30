"""Unique and verified Length/TM/GC/formamid/na/off-targets/kmers/deltaG."""

import glob
import os
import subprocess
import time

import pandas as pd


def pytest_sessionstart(session):
    """Delete old files in case of rerun."""
    files = glob.glob("./*/sacCer3*.ebwt")
    files.extend(glob.glob("./*/sacCer3*.jf"))
    files.extend(glob.glob("./*/Saccharomyces*"))
    files.extend(glob.glob("./*/YKL185W_*.jf"))
    files.extend(glob.glob("./*/renilla_*.jf"))
    [os.remove(f) for f in files]  # type: ignore


def verify_csv(fname: str, *params):
    """Basic inspection if CSV file has the right columns and values fall within ranges."""
    df = pd.read_csv(fname)
    assert all(df["length"] == df["sequence"].apply(len))

    # Check length if provided
    if "--min-length" in params and "--max-length" in params:
        idx_min = params.index("--min-length") + 1
        idx_max = params.index("--max-length") + 1
        assert all(df["length"] >= params[idx_min])
        assert all(df["length"] <= params[idx_max])

    # Check TM if provided
    if "--min-tm" in params and "--max-tm" in params:
        idx_min = params.index("--min-tm") + 1
        idx_max = params.index("--max-tm") + 1
        assert all(df["TM"] >= params[idx_min])
        assert all(df["TM"] <= params[idx_max])

    # Check GC if provided
    if "--min-gc" in params and "--max-gc" in params:
        idx_min = params.index("--min-gc") + 1
        idx_max = params.index("--max-gc") + 1
        assert all(df["GC"] >= params[idx_min])
        assert all(df["GC"] <= params[idx_max])


def test_help_only():
    """Help only with and without being explicit."""
    process_1 = subprocess.check_output("efishent").decode()
    process_2 = subprocess.check_output("eFISHent").decode()
    process_3 = subprocess.check_output(["efishent", "--help"]).decode()
    assert process_1 == process_2 == process_3


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
    bowtie_files = glob.glob("./tests/sacCer3*.ebwt")
    assert len(bowtie_files) == 6
    assert os.path.isfile("./tests/sacCer3_15.jf")


def test_passed_exogenous():
    """Exogenous sequence (renilla) as sequence file, unique output dir, all files saved."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--save-intermediates",
        "True",
        "--sequence-file",
        "./tests/renilla.fasta",
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
    tic = time.time()
    process = subprocess.run(args)
    assert time.time() - tic < 60

    assert process.returncode == 0
    csv_files = glob.glob("./renilla_*.csv")
    assert len(csv_files) == 1
    verify_csv(csv_files[0], args)
    [os.remove(f) for f in glob.glob("./renilla_*")]  # type: ignore


def test_downloaded_endogenous_optimal():
    """Gene/organism downloaded downloaded, minus strand, optimal model w/ timeout."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--gene-name",
        "AAD4",
        "--organism-name",
        "Saccharomyces cerevisiae",
        "--is-plus-strand",
        "False",
        "--optimization-method",
        "optimal",
        "--optimization-time-limit",
        "1",
        "--min-length",
        "45",
        "--max-length",
        "45",
        "--min-tm",
        "61",
        "--max-tm",
        "62",
        "--formamide-concentration",
        "20",
        "--debug",
    ]
    tic = time.time()
    process = subprocess.run(args)
    assert time.time() - tic < 60

    assert process.returncode == 0
    files = sorted(glob.glob("./Saccharomyces_cerevisiae_AAD4_*"))
    assert len(files) == 3
    verify_csv(files[0], args)
    [os.remove(f) for f in files]  # type: ignore


def test_counttable_full_output():
    """Ensembl downloaded sequence with count table and max off target FPKM."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/sacCer3.fa",
        "--reference-annotation",
        "./tests/sacCer3.gtf",
        "--ensembl-id",
        "YKL185W",
        "--encode-count-table",
        "./tests/count_table1.tsv",
        "--max-expression-percentage",
        "40",
        "--debug",
    ]
    tic = time.time()
    process = subprocess.run(args)
    assert time.time() - tic < 60

    assert process.returncode == 0
    files = glob.glob("./YKL185W_*")
    assert len(files) == 3
    [os.remove(f) for f in files]  # type: ignore


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
