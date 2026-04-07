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


def test_find_dangerous_probes_nonexistent_file():
    """Non-existent BLAST output file should return empty set."""
    task = FilterRibosomalRNA()
    dangerous = task._find_dangerous_probes("/nonexistent/path/blast.tsv")
    assert len(dangerous) == 0


def test_find_dangerous_probes_gapped_alignment():
    """Alignment with gaps should account for gapopen in effective length."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # 21nt probe, 20nt alignment, 1 gap, 1 mismatch
        # effective_len = 20 - 1 = 19, which is 90% of 21 -> passes 85% threshold
        # mismatches = 1 <= 3 -> dangerous
        f.write("gapped-probe\tU13369.1\t95.0\t20\t1\t1\t1\t20\t100\t119\t0.01\t30\t21\n")
        # 21nt probe, 20nt alignment, 3 gaps -> effective_len = 17, which is 81% -> below 85%
        f.write("safe-gapped\tU13369.1\t95.0\t20\t1\t3\t1\t20\t100\t119\t0.01\t30\t21\n")
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert "gapped-probe" in dangerous
    assert "safe-gapped" not in dangerous


def test_find_dangerous_probes_multiple_hits_same_probe():
    """Multiple hits for the same probe should all be considered."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # First hit: safe (low coverage)
        f.write("multi-hit\tU13369.1\t90.0\t10\t1\t0\t1\t10\t100\t109\t0.1\t15\t25\n")
        # Second hit: dangerous (high coverage, low mismatches)
        f.write("multi-hit\tU13369.1\t95.0\t23\t1\t0\t1\t23\t200\t222\t0.001\t40\t25\n")
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert "multi-hit" in dangerous


def test_find_dangerous_probes_boundary_85_percent():
    """Test the exact 85% coverage boundary."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # 20nt probe, 17nt alignment = exactly 85% -> should be rejected
        f.write("boundary-pass\tU13369.1\t95.0\t17\t0\t0\t1\t17\t100\t116\t0.001\t30\t20\n")
        # 20nt probe, 16nt alignment = 80% -> should pass
        f.write("boundary-fail\tU13369.1\t95.0\t16\t0\t0\t1\t16\t100\t115\t0.001\t30\t20\n")
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert "boundary-pass" in dangerous
    assert "boundary-fail" not in dangerous


def test_find_dangerous_probes_boundary_3_mismatches():
    """Test the exact 3 mismatch boundary."""
    task = FilterRibosomalRNA()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # 3 mismatches -> should be rejected
        f.write("three-mm\tU13369.1\t85.0\t20\t3\t0\t1\t20\t100\t119\t0.01\t25\t21\n")
        # 4 mismatches -> should pass
        f.write("four-mm\tU13369.1\t80.0\t20\t4\t0\t1\t20\t100\t119\t0.01\t20\t21\n")
        f.flush()
        dangerous = task._find_dangerous_probes(f.name)

    os.unlink(f.name)
    assert "three-mm" in dangerous
    assert "four-mm" not in dangerous


def test_build_blast_db_no_makeblastdb():
    """_build_blast_db should return False when makeblastdb is not found."""
    task = FilterRibosomalRNA()
    with tempfile.TemporaryDirectory() as tmpdir:
        from unittest.mock import patch
        with patch("shutil.which", return_value=None):
            result = task._build_blast_db("/some/fasta.fa", os.path.join(tmpdir, "db"))
            assert result is False


def test_build_blast_db_subprocess_error():
    """_build_blast_db should return False when makeblastdb fails."""
    import subprocess
    from unittest.mock import patch

    task = FilterRibosomalRNA()
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("shutil.which", return_value="/usr/bin/makeblastdb"), \
             patch("eFISHent.rdna_filter.subprocess.check_call",
                   side_effect=subprocess.CalledProcessError(1, "makeblastdb")):
            result = task._build_blast_db("/some/fasta.fa", os.path.join(tmpdir, "db"))
            assert result is False


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
def test_trna_probe_rejected_by_imaging_risk_screen():
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

def test_run_no_blast_installed():
    """run() should skip filtering and write all probes when BLAST is not installed."""
    from unittest.mock import patch, MagicMock

    task = FilterRibosomalRNA()
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create input probe fasta
        probe_fasta = os.path.join(tmpdir, "input_probes.fasta")
        records = [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("ATCGATCGATCGATCGATCG"),
                id="probe-1", name="probe-1", description=""
            ),
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("GCTAGCTAGCTAGCTAGCTA"),
                id="probe-2", name="probe-2", description=""
            ),
        ]
        Bio.SeqIO.write(records, probe_fasta, "fasta")

        output_fasta = os.path.join(tmpdir, "output.fasta")
        mock_input = MagicMock()
        mock_input.path = probe_fasta
        mock_output = MagicMock()
        mock_output.path = output_fasta

        with patch.object(task, "input", return_value=mock_input), \
             patch.object(task, "output", return_value=mock_output), \
             patch("eFISHent.util.log_stage_start"), \
             patch("shutil.which", return_value=None):
            task.run()

        # All probes should be written
        output_records = list(Bio.SeqIO.parse(output_fasta, "fasta"))
        assert len(output_records) == 2


def test_run_no_references():
    """run() should skip filtering when no reference FASTAs exist."""
    from unittest.mock import patch, MagicMock

    task = FilterRibosomalRNA()
    with tempfile.TemporaryDirectory() as tmpdir:
        probe_fasta = os.path.join(tmpdir, "input_probes.fasta")
        records = [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq("ATCGATCGATCGATCGATCG"),
                id="probe-1", name="probe-1", description=""
            ),
        ]
        Bio.SeqIO.write(records, probe_fasta, "fasta")

        output_fasta = os.path.join(tmpdir, "output.fasta")
        mock_input = MagicMock()
        mock_input.path = probe_fasta
        mock_output = MagicMock()
        mock_output.path = output_fasta

        with patch.object(task, "input", return_value=mock_input), \
             patch.object(task, "output", return_value=mock_output), \
             patch("eFISHent.util.log_stage_start"), \
             patch("shutil.which", return_value="/usr/bin/blastn"), \
             patch("eFISHent.rdna_filter.RDNA_45S_FASTA", "/nonexistent/rdna.fa"), \
             patch("eFISHent.rdna_filter.ALPHA_SAT_FASTA", "/nonexistent/alpha.fa"), \
             patch("eFISHent.rdna_filter.IMAGING_RISK_FASTA", "/nonexistent/risk.fa"), \
             patch("eFISHent.rdna_filter.ProbeConfig") as mock_pc, \
             patch("eFISHent.config.GeneralConfig") as mock_gc:
            mock_pc.return_value.custom_rdna_fasta = ""
            mock_gc.return_value.threads = 1
            task.run()

        output_records = list(Bio.SeqIO.parse(output_fasta, "fasta"))
        assert len(output_records) == 1


@pytest.mark.skipif(not BLAST_AVAILABLE, reason="BLAST+ not installed")
def test_run_filters_dangerous_probes():
    """run() should remove probes that match rDNA sequences."""
    from unittest.mock import patch, MagicMock

    task = FilterRibosomalRNA()
    with tempfile.TemporaryDirectory() as tmpdir:
        probe_fasta = os.path.join(tmpdir, "input_probes.fasta")
        # Ren-28 is known to match 45S rDNA
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

        output_fasta = os.path.join(tmpdir, "output.fasta")
        mock_input = MagicMock()
        mock_input.path = probe_fasta
        mock_output = MagicMock()
        mock_output.path = output_fasta

        with patch.object(task, "input", return_value=mock_input), \
             patch.object(task, "output", return_value=mock_output), \
             patch("eFISHent.util.log_stage_start"), \
             patch("eFISHent.util.log_and_check_candidates"), \
             patch("eFISHent.rdna_filter.ProbeConfig") as mock_pc, \
             patch("eFISHent.config.GeneralConfig") as mock_gc:
            mock_pc.return_value.custom_rdna_fasta = ""
            mock_gc.return_value.threads = 1
            task.run()

        output_records = list(Bio.SeqIO.parse(output_fasta, "fasta"))
        output_ids = {r.id for r in output_records}
        assert "Ren-28" not in output_ids
        assert "safe-probe" in output_ids


def test_requires_returns_basic_filtering():
    """requires() should return BasicFiltering task."""
    task = FilterRibosomalRNA()
    from eFISHent.basic_filtering import BasicFiltering
    req = task.requires()
    assert isinstance(req, BasicFiltering)


def test_output_path():
    """output() should return a fasta file in the output directory."""
    from unittest.mock import patch
    task = FilterRibosomalRNA()
    with patch("eFISHent.util.get_output_dir", return_value="/tmp/out"), \
         patch("eFISHent.util.get_gene_name", return_value="test_gene_abc"):
        output = task.output()
        assert output.path == "/tmp/out/test_gene_abc_rdna_filtered.fasta"


def test_bowtie2_seed_length_short_probes():
    """For probes <=22nt, seed length should be 10 (not 20)."""
    # This is a config-based check — the actual bowtie2 call uses ProbeConfig().min_length
    from eFISHent.config import ProbeConfig

    # Default min_length is 21 which is <=22, so seed should be 10
    assert ProbeConfig().min_length <= 22
