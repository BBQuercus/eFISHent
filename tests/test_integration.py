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
        "./tests/data/sacCer3.fa",
        "--build-indices",
        "True",
    ]
    subprocess.run(args, check=True)
    bowtie_files = glob.glob("./tests/data/sacCer3*.ebwt")
    assert len(bowtie_files) == 6
    assert os.path.isfile("./tests/data/sacCer3_15.jf")


def test_passed_exogenous():
    """Exogenous sequence (renilla) as sequence file, unique output dir, all files saved."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/data/sacCer3.fa",
        "--save-intermediates",
        "True",
        "--sequence-file",
        "./tests/data/renilla.fasta",
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
        "--sequence-similarity",
        "95",
        "--debug",
    ]
    subprocess.run(args, check=True)

    csv_files = sorted(glob.glob("./renilla_*.csv"))
    assert len(csv_files) == 2  # _counts and output
    verify_csv(csv_files[0], args)
    [os.remove(f) for f in glob.glob("./renilla_*")]  # type: ignore


def test_downloaded_endogenous_optimal():
    """Gene/organism downloaded downloaded, minus strand, optimal model w/ timeout."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/data/sacCer3.fa",
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
    subprocess.run(args, check=True)
    assert time.time() - tic < 60

    files = sorted(glob.glob("./Saccharomyces_cerevisiae_AAD4_*"))
    assert len(files) == 4  # log, csv, fasta-probes, fasta-ensembl
    verify_csv(files[0], args)
    [os.remove(f) for f in files]  # type: ignore


def test_counttable_full_output():
    """Ensembl downloaded sequence with count table and max off target FPKM."""
    args = [
        "efishent",
        "-g",
        "./tests/data/sacCer3.fa",
        "-a",
        "./tests/data/sacCer3.gtf",
        "-id",
        "YKL185W",
        "-ct",
        "./tests/data/count_table1.tsv",
        "-ep",
        "40",
        "--na-concentration",
        "390",
        "--debug",
    ]
    subprocess.run(args, check=True)
    files = glob.glob("./YKL185W_*")
    assert len(files) == 4  # log, csv, fasta-probes, fasta-ensembl
    [os.remove(f) for f in files]  # type: ignore


def test_analyze():
    """Analysis of probes with provided file."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/data/sacCer3.fa",
        "--analyze-probeset",
        "./tests/data/aad4_subset.fasta",
        "--sequence-file",
        "./tests/data/aad4.fasta",
        "--debug",
    ]
    subprocess.run(args, check=True)
    fname = "./aad4_subset_analysis.pdf"
    assert os.path.exists(fname)
    os.remove(fname)


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
