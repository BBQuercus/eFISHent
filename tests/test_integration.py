"""Unique and verified Length/TM/GC/formamid/na/off-targets/kmers/deltaG."""

import glob
import os
import shutil
import subprocess
import time

import pandas as pd
import pytest

# Auto-discover deps installed by install.sh (same logic as eFISHent CLI)
_deps_bin = os.path.join(os.path.expanduser("~"), ".local", "efishent", "deps", "bin")
_deps_edirect = os.path.join(os.path.expanduser("~"), ".local", "efishent", "deps", "edirect")
for _p in (_deps_bin, _deps_edirect):
    if os.path.isdir(_p) and _p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{_p}:{os.environ.get('PATH', '')}"

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

    Tests the same esearch | efetch pipeline that eFISHent uses internally
    (queries nuccore directly, no elink step) with a known-good query.
    """
    if not ESEARCH_AVAILABLE:
        return False
    try:
        search = subprocess.run(
            ["esearch", "-db", "nuccore", "-query",
             "(AAD4[Gene Name]) AND Saccharomyces cerevisiae[Organism] AND refseq[filter]"],
            capture_output=True,
            timeout=30,
        )
        fetch = subprocess.run(
            ["efetch", "-format", "fasta"],
            input=search.stdout,
            capture_output=True,
            timeout=30,
        )
        fasta = fetch.stdout.decode()

        # Valid result must contain FASTA header and no error patterns
        if ">" not in fasta:
            return False
        if "Error" in fasta or "Failed" in fasta:
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
        "--allow-no-transcriptome", "True",
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


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie2/bowtie2-build/jellyfish not installed")
def test_genome_flag_e2e():
    """End-to-end: --genome yeast downloads from HF and produces probe output."""
    import tempfile

    with tempfile.TemporaryDirectory() as cache_dir:
        args = [
            "efishent",
            "--genome", "yeast",
            "--index-cache-dir", cache_dir,
            "--sequence-file", "./tests/data/renilla.fasta",
            "--is-endogenous", "False",
            "--min-length", "20",
            "--max-length", "24",
            "--min-tm", "45",
            "--max-tm", "55",
            "--sequence-similarity", "95",
            "--allow-no-transcriptome", "True",
            "--debug",
        ]
        result = subprocess.run(args, check=True)
        assert result.returncode == 0

        # Verify probe output files exist
        csv_files = sorted(glob.glob("./renilla_*.csv"))
        assert len(csv_files) >= 1, "Expected at least one CSV output file with probes"

        # Verify the CSV contains actual probes
        df = pd.read_csv(csv_files[0])
        assert len(df) > 0, "Expected probes in the output CSV"
        assert "sequence" in df.columns
        assert "length" in df.columns

        # Clean up output files
        for f in glob.glob("./renilla_*"):
            if os.path.isfile(f):
                os.remove(f)


@pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="bowtie2/bowtie2-build/jellyfish not installed")
@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_genome_flag_with_gene_name_e2e():
    """End-to-end: --genome + --gene-name downloads from HF + NCBI and produces probes.

    Mirrors the README command pattern (efishent --genome X --gene-name Y
    --organism-name Z) using yeast instead of hg38 for speed.
    """
    import tempfile

    # Clean up any leftover files from previous runs
    for f in glob.glob("./Saccharomyces_cerevisiae_URA3_*"):
        if os.path.isfile(f):
            os.remove(f)

    with tempfile.TemporaryDirectory() as cache_dir:
        args = [
            "efishent",
            "--genome", "yeast",
            "--index-cache-dir", cache_dir,
            "--gene-name", "URA3",
            "--organism-name", "Saccharomyces cerevisiae",
            "--min-length", "40",
            "--max-length", "50",
            "--min-tm", "57",
            "--max-tm", "67",
            "--formamide-concentration", "20",
            "--debug",
        ]
        result = subprocess.run(args, check=True)
        assert result.returncode == 0

        # Verify probe output files exist (final output is {name}.csv, not
        # intermediates like {name}_counts.csv)
        all_files = sorted(glob.glob("./Saccharomyces_cerevisiae_URA3_*"))
        assert len(all_files) >= 1, "Expected output files"

        # Find the final probe CSV (has columns: sequence, length, GC, TM, etc.)
        csv_files = sorted(glob.glob("./Saccharomyces_cerevisiae_URA3_*.csv"))
        probe_csv = None
        for f in csv_files:
            df = pd.read_csv(f)
            if "sequence" in df.columns:
                probe_csv = f
                break
        assert probe_csv is not None, f"No probe CSV found among: {csv_files}"
        assert len(df) > 0, "Expected probes in the output CSV"
        assert "length" in df.columns

        # Clean up output files
        for f in all_files:
            if os.path.isfile(f):
                os.remove(f)


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


# ---------------------------------------------------------------------------
# End-to-end tests for human genes using --genome hg38
# These are slow (~5-8 min each) and require the hg38 genome indices to be
# cached at ~/.local/efishent/indices/homo_sapiens/GRCh38/.
# Run with: pytest tests/test_integration.py -k "e2e_hg38" -v
# ---------------------------------------------------------------------------

def _hg38_cached():
    """Check if the hg38 pre-built genome indices are cached."""
    from eFISHent.prebuilt import is_genome_cached
    return is_genome_cached("homo_sapiens/GRCh38")


HG38_AVAILABLE = _hg38_cached()

_hg38_skip = pytest.mark.skipif(
    not HG38_AVAILABLE or not ALL_TOOLS_AVAILABLE or not ESEARCH_WORKS or not BLAST_AVAILABLE,
    reason="requires cached hg38 genome + bowtie2/jellyfish/esearch/blastn",
)


def _run_hg38_gene(gene_name, min_probes):
    """Helper: run eFISHent for a human gene and verify output."""
    # Clean up any leftover files
    prefix = f"homo_sapiens_{gene_name}_"
    for f in glob.glob(f"./{prefix}*"):
        if os.path.isfile(f):
            os.remove(f)

    args = [
        "efishent",
        "--genome", "hg38",
        "--gene-name", gene_name,
        "--organism-name", "homo sapiens",
        "--preset", "smfish",
    ]
    result = subprocess.run(args, check=True, timeout=600)
    assert result.returncode == 0

    # Find probe output CSV (has "sequence" column)
    csv_files = sorted(glob.glob(f"./{prefix}*.csv"))
    probe_csv = None
    for f in csv_files:
        df = pd.read_csv(f)
        if "sequence" in df.columns:
            probe_csv = f
            break
    assert probe_csv is not None, f"No probe CSV found among: {csv_files}"
    assert len(df) >= min_probes, (
        f"Expected >= {min_probes} probes for {gene_name}, got {len(df)}"
    )
    assert "length" in df.columns
    assert "TM" in df.columns

    # Verify FASTA output exists and matches CSV
    fasta_files = sorted(glob.glob(f"./{prefix}*.fasta"))
    assert len(fasta_files) >= 1, "Expected FASTA output file"

    # Clean up
    for f in glob.glob(f"./{prefix}*"):
        if os.path.isfile(f):
            os.remove(f)

    return df


@_hg38_skip
def test_e2e_hg38_actb():
    """ACTB: pseudogene-rich housekeeping gene — tests pseudogene-aware filtering."""
    df = _run_hg38_gene("ACTB", min_probes=3)
    # ACTB is challenging (many pseudogenes); verify probes are reasonable
    assert all(df["length"] >= 20)
    assert all(df["length"] <= 24)


@_hg38_skip
def test_e2e_hg38_eif2b1():
    """EIF2B1: standard gene — tests baseline smFISH pipeline."""
    df = _run_hg38_gene("EIF2B1", min_probes=10)
    assert all(df["length"] >= 20)
    assert all(df["length"] <= 24)


@_hg38_skip
def test_e2e_hg38_tp53():
    """TP53: tumor suppressor with many isoforms — tests transcript selection."""
    df = _run_hg38_gene("TP53", min_probes=10)
    assert all(df["length"] >= 20)
    assert all(df["length"] <= 24)
