"""Unique and verified Length/TM/GC/formamid/na/off-targets/kmers/deltaG."""

import glob
import os
import shutil
import subprocess
import time

import pandas as pd
import pytest

# Check for required external dependencies
BOWTIE_AVAILABLE = shutil.which("bowtie") is not None
BOWTIE_BUILD_AVAILABLE = shutil.which("bowtie-build") is not None
JELLYFISH_AVAILABLE = shutil.which("jellyfish") is not None
BOWTIE2_AVAILABLE = shutil.which("bowtie2") is not None
BOWTIE2_BUILD_AVAILABLE = shutil.which("bowtie2-build") is not None
BLAST_AVAILABLE = shutil.which("blastn") is not None
ESEARCH_AVAILABLE = shutil.which("esearch") is not None
GLPK_AVAILABLE = shutil.which("glpsol") is not None


def _esearch_works():
    """Check if esearch actually works (not just installed).

    Tests the full pipeline (esearch | elink | efetch) with a known-good query
    to verify it can actually retrieve FASTA data.
    """
    if not ESEARCH_AVAILABLE:
        return False
    try:
        # Test full pipeline with a simple, reliable query
        search = subprocess.run(
            ["esearch", "-db", "gene", "-query", "human insulin"],
            capture_output=True,
            timeout=15,
        )
        link = subprocess.run(
            ["elink", "-db", "gene", "-target", "nuccore", "-name", "gene_nuccore_refseqrna"],
            input=search.stdout,
            capture_output=True,
            timeout=15,
        )
        fetch = subprocess.run(
            ["efetch", "-format", "fasta"],
            input=link.stdout,
            capture_output=True,
            timeout=15,
        )
        fasta = fetch.stdout.decode()

        # Valid result must contain FASTA header and no error patterns
        if ">" not in fasta:
            return False
        fasta_nospace = fasta.replace(" ", "")
        if "Error" in fasta or "Failed" in fasta:
            return False
        if "Error" in fasta_nospace or "Failed" in fasta_nospace:
            return False
        return True
    except Exception:
        return False


ESEARCH_WORKS = _esearch_works()

# Integration tests require all core tools
ALL_TOOLS_AVAILABLE = all([
    BOWTIE2_AVAILABLE,
    BOWTIE2_BUILD_AVAILABLE,
    JELLYFISH_AVAILABLE,
])

LEGACY_TOOLS_AVAILABLE = all([
    BOWTIE_AVAILABLE,
    BOWTIE_BUILD_AVAILABLE,
    JELLYFISH_AVAILABLE,
])


def pytest_sessionstart(session):
    """Delete old files in case of rerun."""
    files = glob.glob("./*/sacCer3*.ebwt")
    files.extend(glob.glob("./*/sacCer3*.bt2"))
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


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
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


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
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
        "--allow-no-transcriptome", "True",
        "--debug",
    ]
    subprocess.run(args, check=True)

    csv_files = sorted(glob.glob("./renilla_*.csv"))
    assert len(csv_files) == 2  # _counts and output
    verify_csv(csv_files[0], args)
    [os.remove(f) for f in glob.glob("./renilla_*") if os.path.isfile(f)]  # type: ignore


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
@pytest.mark.skipif(not GLPK_AVAILABLE, reason="GLPK solver (glpsol) not installed")
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


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
@pytest.mark.skip(reason="gtfparse uses deprecated polars API (toggle_string_cache)")
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


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
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


@pytest.mark.skipif(not LEGACY_TOOLS_AVAILABLE, reason="bowtie/bowtie-build/jellyfish not installed")
def test_backward_compat_bowtie1():
    """Running with --aligner bowtie produces output (legacy mode)."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/data/sacCer3.fa",
        "--sequence-file",
        "./tests/data/renilla.fasta",
        "--is-endogenous",
        "False",
        "--aligner",
        "bowtie",
        "--min-length",
        "20",
        "--max-length",
        "24",
        "--min-tm",
        "45",
        "--max-tm",
        "55",
        "--debug",
    ]
    subprocess.run(args, check=True)
    csv_files = sorted(glob.glob("./renilla_*.csv"))
    assert len(csv_files) >= 1
    [os.remove(f) for f in glob.glob("./renilla_*")]  # type: ignore


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie2/bowtie2-build/jellyfish not installed")
def test_low_complexity_filters_active():
    """Pipeline with --max-homopolymer-length 4 should complete."""
    args = [
        "efishent",
        "--reference-genome",
        "./tests/data/sacCer3.fa",
        "--sequence-file",
        "./tests/data/renilla.fasta",
        "--is-endogenous",
        "False",
        "--max-homopolymer-length",
        "4",
        "--min-length",
        "20",
        "--max-length",
        "24",
        "--min-tm",
        "45",
        "--max-tm",
        "55",
        "--allow-no-transcriptome", "True",
        "--debug",
    ]
    subprocess.run(args, check=True)
    csv_files = sorted(glob.glob("./renilla_*.csv"))
    assert len(csv_files) >= 1
    [os.remove(f) for f in glob.glob("./renilla_*") if os.path.isfile(f)]  # type: ignore


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
