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
    with pytest.raises(LookupError):
        task_download.fetch_entrez("asldkjfalsdkfalsdbfalsdif")


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_entrez_gene_organism(task_download):
    query = "hr-38 [GENE] drosophila melanogaster [ORGN]"
    assert isinstance(task_download.fetch_entrez(query), str)


@pytest.mark.skipif(not ESEARCH_WORKS, reason="NCBI E-utilities (esearch) not working")
def test_download_entrez_ensembl(task_download):
    vimentin = task_download.fetch_entrez("ENSG00000026025")
    assert isinstance(vimentin, str)
    assert vimentin.startswith(">NM_003380.5 Homo sapiens vimentin (VIM)")


@pytest.mark.parametrize(
    "ensembl,gene,organism,query",
    [
        ("ENSG00000026025", "", "", "(ENSG00000026025)"),
        (
            "ENSG00000128272",
            "gene",
            "homo sapiens",
            "(ENSG00000128272) AND homo sapiens[Organism]",
        ),
    ],
)
def test_get_ensembl_query(ensembl, gene, organism, query, task_download):
    assert task_download.get_entrez_query(ensembl, gene, organism) == query


@pytest.mark.parametrize(
    "gene,organism", [("hr-38", "drosophila melanogaster"), ("vim", "homo sapiens")]
)
def test_get_ensembl_gene_organism(gene, organism, task_download):
    assert (
        task_download.get_entrez_query("", gene, organism)
        == f"({gene}[Gene Name]) AND {organism}[Organism]"
    )


def test_get_ensembl_error(task_download):
    with pytest.raises(ValueError):
        task_download.get_entrez_query("", "", "")


def test_get_sequence_single(task_prepare, sequence):
    output = task_prepare.select_sequence([sequence])
    assert output.seq == sequence.seq
    assert output.id == sequence.id


def test_get_sequence_multiple(task_prepare, sequence, caplog):
    sequences = [sequence for _ in range(10)]

    output = task_prepare.select_sequence(sequences)

    for record in caplog.records:
        assert record.levelname == "WARNING"

    assert output.seq == sequence.seq
    assert output.id == sequence.id


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


class TestEntrezSubprocessErrors:
    """Tests for proper error handling in entrez subprocess calls."""

    def test_elink_failure_raises_error(self, task_download):
        """Test that elink failure raises CalledProcessError.

        This test verifies that adding check=True to elink will properly
        propagate subprocess failures instead of silently continuing.
        """

        def mock_subprocess_run(*args, **kwargs):
            cmd = args[0]
            mock_result = MagicMock()

            if cmd[0] == "esearch":
                mock_result.returncode = 0
                mock_result.stdout = b"<ENTREZ_DIRECT>valid</ENTREZ_DIRECT>"
                return mock_result
            elif cmd[0] == "elink":
                # Simulate elink failure
                if kwargs.get("check", False):
                    raise subprocess.CalledProcessError(1, cmd)
                mock_result.returncode = 1
                mock_result.stdout = b""
                return mock_result
            elif cmd[0] == "efetch":
                mock_result.returncode = 0
                mock_result.stdout = b""
                return mock_result
            return mock_result

        with patch(
            "eFISHent.prepare_sequence.subprocess.run", side_effect=mock_subprocess_run
        ):
            with pytest.raises(subprocess.CalledProcessError):
                task_download.fetch_entrez("test_query")

    def test_efetch_failure_raises_error(self, task_download):
        """Test that efetch failure raises CalledProcessError.

        This test verifies that adding check=True to efetch will properly
        propagate subprocess failures instead of silently continuing.
        """

        def mock_subprocess_run(*args, **kwargs):
            cmd = args[0]
            mock_result = MagicMock()

            if cmd[0] == "esearch":
                mock_result.returncode = 0
                mock_result.stdout = b"<ENTREZ_DIRECT>valid</ENTREZ_DIRECT>"
                return mock_result
            elif cmd[0] == "elink":
                mock_result.returncode = 0
                mock_result.stdout = b"<ENTREZ_DIRECT>linked</ENTREZ_DIRECT>"
                return mock_result
            elif cmd[0] == "efetch":
                # Simulate efetch failure
                if kwargs.get("check", False):
                    raise subprocess.CalledProcessError(1, cmd)
                mock_result.returncode = 1
                mock_result.stdout = b""
                return mock_result
            return mock_result

        with patch(
            "eFISHent.prepare_sequence.subprocess.run", side_effect=mock_subprocess_run
        ):
            with pytest.raises(subprocess.CalledProcessError):
                task_download.fetch_entrez("test_query")
