import shutil
import subprocess
from unittest.mock import MagicMock, patch

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.prepare_sequence import DownloadEntrezGeneSequence
from eFISHent.prepare_sequence import PrepareSequence

ESEARCH_AVAILABLE = shutil.which("esearch") is not None


def _esearch_works():
    """Check if esearch actually works (not just installed).

    Tests a simple nuccore query to verify it can retrieve FASTA data.
    """
    if not ESEARCH_AVAILABLE:
        return False
    try:
        search = subprocess.run(
            [
                "esearch", "-db", "nuccore", "-query",
                "(actb[Gene Name]) AND mus musculus[Organism]"
                " AND refseq[filter] AND (biomol_mrna[prop] OR biomol_rna[prop])",
            ],
            capture_output=True,
            timeout=15,
        )
        fetch = subprocess.run(
            ["efetch", "-format", "fasta"],
            input=search.stdout,
            capture_output=True,
            timeout=15,
        )
        fasta = fetch.stdout.decode()
        return ">" in fasta and "NM_" in fasta
    except Exception:
        return False


ESEARCH_WORKS = _esearch_works()


@pytest.fixture
def task_download():
    return DownloadEntrezGeneSequence()


@pytest.fixture
def task_prepare():
    return PrepareSequence()


@pytest.fixture
def sequence():
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="my_id")


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_entrez_error(task_download):
    with pytest.raises((LookupError, subprocess.CalledProcessError)):
        task_download.fetch_entrez("", "asldkjfalsdkfalsdbf", "nonexistent organism")


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_entrez_gene_organism(task_download):
    fasta = task_download.fetch_entrez("", "actb", "mus musculus")
    assert isinstance(fasta, str)
    assert ">" in fasta
    assert "NM_" in fasta


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_entrez_ensembl(task_download):
    vimentin = task_download.fetch_entrez("ENSG00000026025", "", "")
    assert isinstance(vimentin, str)
    assert "NM_" in vimentin
    assert "vimentin" in vimentin.lower() or "VIM" in vimentin


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_returns_only_refseq(task_download):
    """Verify that downloads return only RefSeq transcripts, not genomic sequences."""
    fasta = task_download.fetch_entrez("", "gapdh", "homo sapiens")
    # All records should be NM_ or NR_ (RefSeq transcripts)
    for line in fasta.split("\n"):
        if line.startswith(">"):
            accession = line.split()[0][1:]
            assert accession.startswith(("NM_", "NR_")), (
                f"Non-RefSeq record in download: {accession}"
            )


def test_get_nuccore_query_gene_organism(task_download):
    query = task_download.get_nuccore_query("", "vim", "homo sapiens")
    assert "(vim[Gene Name])" in query
    assert "homo sapiens[Organism]" in query
    assert "refseq[filter]" in query
    assert "biomol_mrna[prop]" in query


def test_get_nuccore_query_error(task_download):
    with pytest.raises(ValueError):
        task_download.get_nuccore_query("", "", "")


def test_get_sequence_single(task_prepare, sequence):
    output = task_prepare.select_sequence([sequence])
    assert output.seq == sequence.seq
    assert output.id == sequence.id


def test_get_sequence_multiple_selects_longest(task_prepare):
    short = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCG"), id="short")
    long = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCGATCG"), id="long")
    medium = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="medium")

    output = task_prepare.select_sequence([short, long, medium])

    assert output.id == "long"
    assert len(output.seq) == 12


def test_get_sequence_error(task_prepare):
    with pytest.raises(ValueError):
        assert task_prepare.select_sequence([])


def test_get_strand_plus(task_prepare, sequence):
    assert task_prepare.select_strand(sequence, is_plus_strand=True).seq == Bio.Seq.Seq(
        "CGATCGAT"
    )


def test_get_strand_minus(task_prepare, sequence):
    assert task_prepare.select_strand(
        sequence, is_plus_strand=False
    ).seq == Bio.Seq.Seq("ATCGATCG")


def test_properly_formatted_fasta_single(task_prepare, tmpdir):
    fasta_file = tmpdir.join("proper_fasta.fasta")
    fasta_file.write(">Sequence1\nATCGATCGANNN\nAGTGGGGATGTAAA")
    assert task_prepare.is_fasta_formatted(fasta_file) is True


def test_properly_formatted_fasta_multiple(task_prepare, tmpdir):
    fasta_file = tmpdir.join("proper_fasta.fasta")
    fasta_file.write(">Sequence1\nATCGATCGANNN\n>Sequence2\nAGCTAGCTAGCT")
    assert task_prepare.is_fasta_formatted(fasta_file) is True


def test_improperly_formatted_fasta(task_prepare, tmpdir):
    fasta_file = tmpdir.join("improper_fasta.fasta")
    fasta_file.write(">Sequence1\nATCGATCGATCG\nSequence2\nAGCTAGCTAGCT")
    assert task_prepare.is_fasta_formatted(fasta_file) is False


def test_wrong_characters_in_fasta(task_prepare, tmpdir):
    fasta_file = tmpdir.join("improper_fasta.fasta")
    fasta_file.write(">Sequence1\nATCGANNXYZ")
    assert task_prepare.is_fasta_formatted(fasta_file) is False


def test_non_existent_file(task_prepare, tmpdir):
    non_existent_file = str(tmpdir.join("non_existent.fasta"))
    assert task_prepare.is_fasta_formatted(non_existent_file) is False


class TestRegionSelection:
    """Tests for target region selection (exon/intron/both)."""

    def _make_seq(self, seq_str, seq_id="test_gene"):
        return Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(seq_str), id=seq_id, name=seq_id, description=""
        )

    def test_exon_mode_returns_unchanged(self):
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCG" * 10)
        result = task._apply_region_selection(seq, "exon", "TP53")
        assert str(result.seq) == str(seq.seq)

    def test_intron_mode_fallback_without_annotation(self):
        """Intron mode without annotation should return unchanged."""
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCG" * 10)
        result = task._apply_region_selection(seq, "intron", "TP53")
        assert str(result.seq) == str(seq.seq)

    def test_both_mode_fallback_without_annotation(self):
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCG" * 10)
        result = task._apply_region_selection(seq, "both", "TP53")
        assert str(result.seq) == str(seq.seq)


class TestEntrezSubprocessErrors:
    """Tests for proper error handling in entrez subprocess calls."""

    def test_esearch_failure_raises_error(self, task_download):
        """Test that esearch failure raises CalledProcessError."""

        def mock_subprocess_run(*args, **kwargs):
            cmd = args[0]
            if cmd[0] == "esearch":
                raise subprocess.CalledProcessError(1, cmd)
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = b""
            return mock_result

        with patch(
            "eFISHent.prepare_sequence.subprocess.run", side_effect=mock_subprocess_run
        ):
            with pytest.raises(subprocess.CalledProcessError):
                task_download.fetch_entrez("", "test_gene", "test_organism")

    def test_efetch_failure_raises_error(self, task_download):
        """Test that efetch failure raises CalledProcessError."""

        def mock_subprocess_run(*args, **kwargs):
            cmd = args[0]
            mock_result = MagicMock()

            if cmd[0] == "esearch":
                mock_result.returncode = 0
                mock_result.stdout = b"<ENTREZ_DIRECT>valid</ENTREZ_DIRECT>"
                return mock_result
            elif cmd[0] == "efetch":
                raise subprocess.CalledProcessError(1, cmd)
            return mock_result

        with patch(
            "eFISHent.prepare_sequence.subprocess.run", side_effect=mock_subprocess_run
        ):
            with pytest.raises(subprocess.CalledProcessError):
                task_download.fetch_entrez("", "test_gene", "test_organism")

    def test_empty_result_raises_lookup_error(self, task_download):
        """Test that empty FASTA results raise LookupError."""

        def mock_subprocess_run(*args, **kwargs):
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = b"No items found."
            return mock_result

        with patch(
            "eFISHent.prepare_sequence.subprocess.run", side_effect=mock_subprocess_run
        ):
            with pytest.raises(LookupError):
                task_download.fetch_entrez("", "nonexistent", "nonexistent")
