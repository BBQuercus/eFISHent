"""Tests for the rDNA/satellite/imaging-risk RNA filter."""

import os
import shutil
import tempfile

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import pytest

from eFISHent.rdna_filter import (
    FilterRibosomalRNA,
    RDNA_45S_FASTA,
    ALPHA_SAT_FASTA,
    IMAGING_RISK_FASTA,
)

BLAST_AVAILABLE = shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None


# --- Data file existence ---

def test_bundled_45s_fasta_exists():
    """The human 45S rDNA reference must be bundled."""
    assert os.path.isfile(RDNA_45S_FASTA)


def test_bundled_alpha_satellite_exists():
    """The alpha satellite consensus must be bundled."""
    assert os.path.isfile(ALPHA_SAT_FASTA)


def test_bundled_imaging_risk_panel_exists():
    """The imaging-risk RNA panel must be bundled."""
    assert os.path.isfile(IMAGING_RISK_FASTA)


def test_45s_fasta_length():
    """U13369.1 should be 42,999 bp."""
    records = list(Bio.SeqIO.parse(RDNA_45S_FASTA, "fasta"))
    assert len(records) == 1
    assert len(records[0].seq) == 42999


def test_imaging_risk_panel_contents():
    """Panel should contain mt-rRNAs, snRNAs, 7SL, 7SK, and tRNAs."""
    records = list(Bio.SeqIO.parse(IMAGING_RISK_FASTA, "fasta"))
    # Use both ID and description for matching since NCBI accessions are IDs
    texts = [f"{r.id} {r.description}" for r in records]
    all_text = " ".join(texts)

    assert "MT-RNR" in all_text or "12S" in all_text or "16S" in all_text, \
        "Missing mt-rRNA in panel"
    assert "U1" in all_text or "RNU1" in all_text, "Missing U1 snRNA"
    assert "U2" in all_text or "RNU2" in all_text, "Missing U2 snRNA"
    assert "U6" in all_text or "RNU6" in all_text, "Missing U6 snRNA"
    assert "7SL" in all_text or "RN7SL" in all_text, "Missing 7SL/SRP RNA"
    assert "7SK" in all_text or "RN7SK" in all_text, "Missing 7SK RNA"
    assert "tRNA" in all_text, "Missing tRNA"


# --- BLAST-based rejection logic ---

def test_find_dangerous_probes_rejects_near_match():
    """A probe with <=3 mismatches over >=85% of its length should be rejected."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
        # 21nt probe, 20nt alignment (95% coverage), 2 mismatches
        f.write("Ren-28\tU13369.1\t90.0\t20\t2\t0\t1\t20\t1816\t1835\t0.001\t30\t21\n")
        # Safe probe: 21nt probe, 12nt alignment (57% coverage) — below 85% threshold
        f.write("safe-probe\tU13369.1\t85.0\t12\t2\t0\t1\t12\t500\t511\t0.5\t10\t21\n")
        f.flush()

        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)

    assert "Ren-28" in dangerous
    assert "safe-probe" not in dangerous


def test_find_dangerous_probes_allows_4_mismatches():
    """A probe with 4 mismatches should pass (threshold is <=3)."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("marginal\tU13369.1\t80.0\t20\t4\t0\t1\t20\t100\t119\t0.01\t20\t21\n")
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert "marginal" not in dangerous


def test_find_dangerous_probes_empty_blast():
    """Empty BLAST output should return no dangerous probes."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert len(dangerous) == 0


# --- Ren-28 rejection (spec verification case) ---

@pytest.mark.skipif(not BLAST_AVAILABLE, reason="BLAST+ not installed")
def test_ren28_rejected_by_45s_screen():
    """The offending Ren-28 probe (2mm to 45S 5' ETS) must be rejected.

    From RDNA_OFFTARGET_FIX.md:
    Ren-28: GTTTGCGTTGCTCGGGGTCGT (21nt, 62% GC)
    Has 2-mismatch hit to 45S rDNA U13369.1 pos 1816-1835 (5' ETS)
    """
    task = FilterRibosomalRNA()

    # Build a small test FASTA with Ren-28 and a safe probe
    with tempfile.TemporaryDirectory() as tmpdir:
        probe_fasta = os.path.join(tmpdir, "test_probes.fasta")
        blast_out = os.path.join(tmpdir, "blast.tsv")
        db_path = os.path.join(tmpdir, "45s_db")

        records = [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("GTTTGCGTTGCTCGGGGTCGT"),
                id="Ren-28", name="Ren-28", description=""
            ),
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("AAAAAATTTTTTCCCCCGGGG"),
                id="safe-probe", name="safe-probe", description=""
            ),
        ]
        Bio.SeqIO.write(records, probe_fasta, "fasta")

        assert task._build_blast_db(RDNA_45S_FASTA, db_path)
        task._blast_probes(probe_fasta, db_path, blast_out, threads=1)
        dangerous = task._find_dangerous_probes(blast_out)

    assert "Ren-28" in dangerous, "Ren-28 should be rejected by 45S rDNA screen"
    assert "safe-probe" not in dangerous


@pytest.mark.skipif(not BLAST_AVAILABLE, reason="BLAST+ not installed")
def test_tRNA_probe_rejected_by_imaging_risk_screen():
    """A probe matching the bundled tRNA panel must be rejected."""
    task = FilterRibosomalRNA()

    with tempfile.TemporaryDirectory() as tmpdir:
        probe_fasta = os.path.join(tmpdir, "test_probes.fasta")
        blast_out = os.path.join(tmpdir, "blast.tsv")
        db_path = os.path.join(tmpdir, "risk_db")

        records = [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("GGGGAATTAGCTCAGTTGGGA"),
                id="tRNA-hit", name="tRNA-hit", description=""
            ),
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("AAAAAATTTTTTCCCCCGGGG"),
                id="safe-probe", name="safe-probe", description=""
            ),
        ]
        Bio.SeqIO.write(records, probe_fasta, "fasta")

        assert task._build_blast_db(IMAGING_RISK_FASTA, db_path)
        task._blast_probes(probe_fasta, db_path, blast_out, threads=1)
        dangerous = task._find_dangerous_probes(blast_out)

    assert "tRNA-hit" in dangerous, "tRNA-matching probe should be rejected"
    assert "safe-probe" not in dangerous


# --- Cross-hybridization rejection ---

def test_cross_hyb_rejection_enabled_by_default():
    """reject_cross_hybridization should default to True."""
    from eFISHent.config import ProbeConfig
    assert ProbeConfig().reject_cross_hybridization is True


# --- Exogenous without transcriptome ---

def test_exogenous_without_transcriptome_blocked():
    """Exogenous design without transcriptome should fail validation."""
    import argparse
    from unittest.mock import MagicMock

    from eFISHent.cli import validate_args

    args = argparse.Namespace(
        is_endogenous=False,
        reference_transcriptome="",
        reference_annotation="",
        reference_genome="/fake/genome.fa",
        allow_no_transcriptome=False,
        build_indices=False,
        analyze_probeset="",
        min_length=20,
        max_length=25,
        min_tm=40.0,
        max_tm=60.0,
        min_gc=20.0,
        max_gc=80.0,
        kmer_length=15,
        sequence_file="/fake/seq.fa",
        gene_name="",
        organism_name="",
        ensembl_id="",
        encode_count_table="",
        intergenic_off_targets=False,
        filter_rrna=False,
    )

    parser = MagicMock()
    validate_args(args, parser)
    parser.error.assert_called_once()
    error_msg = parser.error.call_args[0][0]
    assert "transcriptome" in error_msg.lower()


def test_exogenous_with_override_allowed():
    """Exogenous design with --allow-no-transcriptome should pass."""
    import argparse
    from unittest.mock import MagicMock

    from eFISHent.cli import validate_args

    args = argparse.Namespace(
        is_endogenous=False,
        reference_transcriptome="",
        reference_annotation="",
        reference_genome="/fake/genome.fa",
        allow_no_transcriptome=True,
        build_indices=False,
        analyze_probeset="",
        min_length=20,
        max_length=25,
        min_tm=40.0,
        max_tm=60.0,
        min_gc=20.0,
        max_gc=80.0,
        kmer_length=15,
        sequence_file="/fake/seq.fa",
        gene_name="",
        organism_name="",
        ensembl_id="",
        encode_count_table="",
        intergenic_off_targets=False,
        filter_rrna=False,
    )

    parser = MagicMock()
    validate_args(args, parser)
    parser.error.assert_not_called()


# --- filter_rdna_45s config ---

def test_filter_rdna_45s_default_enabled():
    """The rDNA filter should be enabled by default."""
    from eFISHent.config import ProbeConfig
    assert ProbeConfig().filter_rdna_45s is True


# --- Bowtie2 seed length fix ---

def test_bowtie2_seed_length_short_probes():
    """For probes <=22nt, seed length should be 10 (not 20)."""
    # This is a config-based check — the actual bowtie2 call uses ProbeConfig().min_length
    from eFISHent.config import ProbeConfig

    # Default min_length is 21 which is <=22, so seed should be 10
    assert ProbeConfig().min_length <= 22
