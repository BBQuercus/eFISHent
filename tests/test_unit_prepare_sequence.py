import shutil
import subprocess
from unittest.mock import MagicMock, patch

import Bio.Seq
import Bio.SeqRecord
import pandas as pd
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


class TestIsFastaFormattedEdgeCases:
    """Additional edge cases for is_fasta_formatted."""

    def test_empty_file(self, tmpdir):
        """Empty file should be considered valid (no invalid content)."""
        task = PrepareSequence()
        fasta_file = tmpdir.join("empty.fasta")
        fasta_file.write("")
        # An empty file has no content to invalidate, returns True
        assert task.is_fasta_formatted(str(fasta_file)) is True

    def test_only_header(self, tmpdir):
        """File with only a header should be valid."""
        task = PrepareSequence()
        fasta_file = tmpdir.join("header_only.fasta")
        fasta_file.write(">Sequence1")
        assert task.is_fasta_formatted(str(fasta_file)) is True

    def test_blank_lines_ignored(self, tmpdir):
        """Blank lines should be ignored."""
        task = PrepareSequence()
        fasta_file = tmpdir.join("blanks.fasta")
        fasta_file.write(">Seq\n\nATCG\n\nGGGG\n")
        assert task.is_fasta_formatted(str(fasta_file)) is True

    def test_lowercase_bases(self, tmpdir):
        """Lowercase bases should be accepted."""
        task = PrepareSequence()
        fasta_file = tmpdir.join("lower.fasta")
        fasta_file.write(">Seq\natcgatcg\n")
        assert task.is_fasta_formatted(str(fasta_file)) is True


class TestSelectSequenceEdgeCases:
    """Additional tests for select_sequence."""

    def test_single_sequence_returned_directly(self):
        """Single sequence should be returned without modification."""
        task = PrepareSequence()
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCG"), id="only")
        result = task.select_sequence([seq])
        assert result is seq

    def test_two_sequences_selects_longest(self):
        """With two sequences, the longer one should be selected."""
        task = PrepareSequence()
        short = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("AT"), id="short")
        long_seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCG"), id="long")
        result = task.select_sequence([short, long_seq])
        assert result.id == "long"

    def test_multiple_sequences_filter_to_gene_name_before_selecting_longest(self):
        """Related-gene transcripts should be excluded before longest-isoform selection."""
        task = PrepareSequence()
        actb_short = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATCG"),
            id="actb-short",
            description="NM_001 ACTB transcript variant 1",
        )
        actb_long = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATCGATCG"),
            id="actb-long",
            description="NM_002 actb transcript variant 2",
        )
        unrelated_longer = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATCGATCGATCG"),
            id="potef",
            description="NM_003 POTEF transcript variant 1",
        )

        with patch("eFISHent.prepare_sequence.SequenceConfig") as mock_cfg:
            mock_cfg.return_value.gene_name = "ACTB"
            result = task.select_sequence([actb_short, actb_long, unrelated_longer])

        assert result.id == "actb-long"

    def test_multiple_sequences_without_gene_match_still_selects_longest(self):
        """If no descriptions match the gene name, fall back to longest sequence."""
        task = PrepareSequence()
        shorter = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATCG"),
            id="tx1",
            description="NM_001 hypothetical transcript",
        )
        longer = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq("ATCGATCG"),
            id="tx2",
            description="NM_002 hypothetical transcript",
        )

        with patch("eFISHent.prepare_sequence.SequenceConfig") as mock_cfg:
            mock_cfg.return_value.gene_name = "ACTB"
            result = task.select_sequence([shorter, longer])

        assert result.id == "tx2"


class TestGetNuccoreQuery:
    """Tests for nuccore query building."""

    def test_ensembl_id_query_needs_esearch(self, task_download):
        """Ensembl ID query should call esearch for gene ID resolution."""
        with patch("eFISHent.prepare_sequence.subprocess.run") as mock_run:
            mock_search = MagicMock()
            mock_search.stdout = b"<ENTREZ_DIRECT>xml</ENTREZ_DIRECT>"
            mock_fetch = MagicMock()
            mock_fetch.stdout = b"12345"
            mock_run.side_effect = [mock_search, mock_fetch]

            query = task_download.get_nuccore_query("ENSG00000026025", "", "")
            assert "12345[Gene ID]" in query
            assert "refseq[filter]" in query

    def test_ensembl_with_organism(self, task_download):
        """Ensembl ID with organism should include organism in esearch query."""
        with patch("eFISHent.prepare_sequence.subprocess.run") as mock_run:
            mock_search = MagicMock()
            mock_search.stdout = b"<ENTREZ_DIRECT>xml</ENTREZ_DIRECT>"
            mock_fetch = MagicMock()
            mock_fetch.stdout = b"12345"
            mock_run.side_effect = [mock_search, mock_fetch]

            task_download.get_nuccore_query("ENSG00000026025", "", "homo sapiens")
            # Check esearch was called with organism
            esearch_cmd = mock_run.call_args_list[0][0][0]
            query_str = esearch_cmd[esearch_cmd.index("-query") + 1]
            assert "homo sapiens" in query_str

    def test_ensembl_no_gene_found_raises(self, task_download):
        """Empty gene ID result should raise LookupError."""
        with patch("eFISHent.prepare_sequence.subprocess.run") as mock_run:
            mock_search = MagicMock()
            mock_search.stdout = b"<ENTREZ_DIRECT>xml</ENTREZ_DIRECT>"
            mock_fetch = MagicMock()
            mock_fetch.stdout = b""
            mock_run.side_effect = [mock_search, mock_fetch]

            with pytest.raises(LookupError, match="No NCBI Gene found"):
                task_download.get_nuccore_query("ENSG_INVALID", "", "")

    def test_ensembl_multiple_gene_ids_uses_first(self, task_download):
        """Multiple gene IDs should use the first one."""
        with patch("eFISHent.prepare_sequence.subprocess.run") as mock_run:
            mock_search = MagicMock()
            mock_search.stdout = b"<ENTREZ_DIRECT>xml</ENTREZ_DIRECT>"
            mock_fetch = MagicMock()
            mock_fetch.stdout = b"12345\n67890\n"
            mock_run.side_effect = [mock_search, mock_fetch]

            query = task_download.get_nuccore_query("ENSG00000026025", "", "")
            assert "12345[Gene ID]" in query


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

    def test_intron_mode_with_annotation_extracts_introns(self):
        """Intron mode with annotation and exons should extract intronic regions."""
        task = PrepareSequence()
        # Sequence: EXON1(0-10) INTRON(10-20) EXON2(20-30)
        seq_str = "AAAAAAAAAA" + "CCCCCCCCCC" + "GGGGGGGGGG"
        seq = self._make_seq(seq_str)
        exons = [(0, 10), (20, 30)]

        with patch.object(task, "_extract_exon_boundaries", return_value=exons), \
             patch("eFISHent.prepare_sequence.GeneralConfig") as mock_gc:
            mock_gc.return_value.reference_annotation = "/fake/annotation.gtf"
            mock_gc.return_value.reference_genome = "/fake/genome.fa"
            result = task._apply_region_selection(seq, "intron", "TP53")
            assert str(result.seq) == "CCCCCCCCCC"
            assert "introns" in result.id

    def test_intron_mode_single_exon_fallback(self):
        """Intron mode for single-exon gene should return unchanged."""
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCGATCGATCG")
        exons = [(0, 16)]  # Single exon covering full length

        with patch.object(task, "_extract_exon_boundaries", return_value=exons), \
             patch("eFISHent.prepare_sequence.GeneralConfig") as mock_gc:
            mock_gc.return_value.reference_annotation = "/fake/annotation.gtf"
            mock_gc.return_value.reference_genome = "/fake/genome.fa"
            result = task._apply_region_selection(seq, "intron", "TP53")
            # Should return unchanged since no intronic regions
            assert str(result.seq) == "ATCGATCGATCGATCG"

    def test_both_mode_with_annotation(self):
        """Both mode should return the full pre-mRNA sequence."""
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCGATCG")
        exons = [(0, 4), (8, 12)]

        with patch.object(task, "_extract_exon_boundaries", return_value=exons), \
             patch("eFISHent.prepare_sequence.GeneralConfig") as mock_gc:
            mock_gc.return_value.reference_annotation = "/fake/annotation.gtf"
            mock_gc.return_value.reference_genome = "/fake/genome.fa"
            result = task._apply_region_selection(seq, "both", "TP53")
            assert str(result.seq) == "ATCGATCGATCG"

    def test_intron_mode_no_exons_found_fallback(self):
        """Intron mode with no exons found should return unchanged."""
        task = PrepareSequence()
        seq = self._make_seq("ATCGATCG" * 10)

        with patch.object(task, "_extract_exon_boundaries", return_value=None), \
             patch("eFISHent.prepare_sequence.GeneralConfig") as mock_gc:
            mock_gc.return_value.reference_annotation = "/fake/annotation.gtf"
            mock_gc.return_value.reference_genome = "/fake/genome.fa"
            result = task._apply_region_selection(seq, "intron", "TP53")
            assert str(result.seq) == str(seq.seq)


class TestExtractExonBoundaries:
    """Tests for _extract_exon_boundaries method."""

    def test_parse_error_returns_none(self):
        """If GTF parsing fails, should return None."""
        task = PrepareSequence()
        # Neither parquet nor raw file exists → returns None
        with patch("os.path.isfile", return_value=False):
            result = task._extract_exon_boundaries("TP53", "/fake/annotation.gtf")
            assert result is None

    def test_no_matching_gene_returns_none(self):
        """If gene not found in GTF, should return None."""
        task = PrepareSequence()
        df = pd.DataFrame({
            "feature": ["exon"],
            "gene_name": ["OTHER_GENE"],
            "gene_id": ["ENSG99999"],
            "start": [100],
            "end": [200],
            "transcript_id": ["ENST00001"],
        })
        with patch("os.path.isfile", return_value=True), \
             patch("pandas.read_parquet", return_value=df):
            result = task._extract_exon_boundaries("TP53", "/fake/annotation.gtf")
            assert result is None

    def test_matching_gene_returns_exons(self):
        """Should return sorted exon boundaries in relative coordinates."""
        task = PrepareSequence()
        df = pd.DataFrame({
            "feature": ["exon", "exon", "gene"],
            "gene_name": ["TP53", "TP53", "TP53"],
            "gene_id": ["ENSG001", "ENSG001", "ENSG001"],
            "start": [1000, 2000, 1000],
            "end": [1100, 2200, 2200],
            "transcript_id": ["ENST001", "ENST001", "ENST001"],
        })
        with patch("os.path.isfile", return_value=True), \
             patch("pandas.read_parquet", return_value=df):
            result = task._extract_exon_boundaries("TP53", "/fake/annotation.gtf")
            assert result is not None
            assert len(result) == 2
            assert result[0] == (0, 100)  # 1000-1000, 1100-1000
            assert result[1] == (1000, 1200)  # 2000-1000, 2200-1000

    def test_gene_id_fallback(self):
        """If gene_name doesn't match, should try gene_id."""
        task = PrepareSequence()
        df = pd.DataFrame({
            "feature": ["exon"],
            "gene_name": ["SomeOtherName"],
            "gene_id": ["TP53_gene"],
            "start": [500],
            "end": [600],
            "transcript_id": ["ENST001"],
        })
        with patch("os.path.isfile", return_value=True), \
             patch("pandas.read_parquet", return_value=df):
            result = task._extract_exon_boundaries("TP53", "/fake/annotation.gtf")
            assert result is not None
            assert len(result) == 1


class TestPrepareSequenceOutput:
    """Tests for PrepareSequence output and requires."""

    def test_output_path(self):
        """output() should return a fasta file in the output directory."""
        task = PrepareSequence()
        with patch("eFISHent.util.get_gene_name", return_value="test_gene_abc"), \
             patch("eFISHent.util.get_output_dir", return_value="/tmp/out"):
            output = task.output()
            assert output.path == "/tmp/out/test_gene_abc_sequence.fasta"

    def test_requires_with_sequence_file(self):
        """With a sequence_file, requires() should return empty dict."""
        task = PrepareSequence()
        with patch("eFISHent.prepare_sequence.SequenceConfig") as mock_sc:
            mock_sc.return_value.sequence_file = "/fake/seq.fasta"
            reqs = task.requires()
            assert reqs == {}

    def test_requires_without_sequence_file(self):
        """Without a sequence_file, requires() should include entrez download."""
        task = PrepareSequence()
        with patch("eFISHent.prepare_sequence.SequenceConfig") as mock_sc:
            mock_sc.return_value.sequence_file = ""
            reqs = task.requires()
            assert "entrez" in reqs


class TestDownloadEntrezOutput:
    """Tests for DownloadEntrezGeneSequence output."""

    def test_output_path(self):
        """output() should return an _entrez.fasta file in the output directory."""
        task = DownloadEntrezGeneSequence()
        with patch("eFISHent.util.get_gene_name", return_value="test_gene_abc"), \
             patch("eFISHent.util.get_output_dir", return_value="/tmp/out"):
            output = task.output()
            assert output.path == "/tmp/out/test_gene_abc_entrez.fasta"


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
