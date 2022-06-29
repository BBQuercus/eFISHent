import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.prepare_sequence import DownloadEntrezGeneSequence
from eFISHent.prepare_sequence import PrepareSequence


@pytest.fixture
def task_download():
    return DownloadEntrezGeneSequence()


@pytest.fixture
def task_prepare():
    return PrepareSequence()


@pytest.fixture
def sequence():
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="my_id")


def test_download_entrez_error(task_download):
    with pytest.raises(LookupError):
        task_download.fetch_entrez("asldkjfalsdkfalsdbfalsdif")


def test_download_entrez_gene_organism(task_download):
    query = "hr-38 [GENE] drosophila melanogaster [ORGN]"
    assert isinstance(task_download.fetch_entrez(query), str)


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
