import configparser
import os
import shutil
import tempfile

import Bio.SeqIO
import luigi
import pytest

from eFISHent.util import create_config_hash
from eFISHent.util import create_data_table
from eFISHent.util import get_gene_name
from eFISHent.util import get_genome_name
from eFISHent.util import get_output_dir
from eFISHent.util import log_and_check_candidates
from eFISHent.util import secure_filename


def test_create_config_hash():
    # Original
    luigi.configuration.add_config_path("./eFISHent/luigi.cfg")
    hash1 = create_config_hash(luigi.configuration.get_config())

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Identical in different place
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        shutil.copy("./eFISHent/luigi.cfg", config_file)
        luigi.configuration.add_config_path(config_file)
        hash2 = create_config_hash(luigi.configuration.get_config())

        # Single modification
        config = configparser.ConfigParser()
        config.read(config_file)
        config.set(
            "GeneralConfig", "threads", str(int(config["GeneralConfig"]["threads"]) + 1)
        )
        with open(config_file, "w") as f:
            config.write(f)
        luigi.configuration.add_config_path(config_file)
        hash3 = create_config_hash(luigi.configuration.get_config())

    assert hash1 == hash2
    assert hash1 != hash3


def test_create_data_table():
    sequences = list(Bio.SeqIO.parse("./tests/data/renilla_data_table.fasta", "fasta"))
    df = create_data_table(sequences)

    for col in ["name", "length", "sequence", "start", "end"]:
        assert col in df.columns

    raw_sequences = [seq.seq for seq in sequences]
    for seq in df["sequence"]:
        assert seq in raw_sequences


def test_gene_name_genefile():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("./tests/data/renilla.fasta")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("")
        organism_name = luigi.Parameter("")

    assert get_gene_name(config=Config).startswith("renilla")
    assert get_gene_name(hashed=False, config=Config) == "renilla"


def test_gene_name_ensembl():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("")
        ensembl_id = luigi.Parameter("ENSG00000026025")
        gene_name = luigi.Parameter("")
        organism_name = luigi.Parameter("")

    assert get_gene_name(config=Config).startswith("ENSG00000026025")


def test_gene_name_genename():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("ACOOLGene123")
        organism_name = luigi.Parameter("Latin name")

    assert get_gene_name(config=Config).startswith("Latin_name_ACOOLGene123_")


def test_gene_name_genefile_genename():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("./tests/data/renilla.fasta")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("ACOOLGene123")
        organism_name = luigi.Parameter("Latin name")

    assert get_gene_name(config=Config).startswith("renilla")


def test_genome_name_notpassed():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_invalid():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("./some_random_file.notfasta")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_nonexistent_with_valid_extension():
    """Test that non-existent file with valid FASTA extension raises ValueError.

    This test exposes the bug where the validation uses OR instead of AND,
    allowing non-existent files with .fa extension to pass validation.
    """
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("./nonexistent_genome.fa")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_exists_with_invalid_extension():
    """Test that existing file with invalid extension raises ValueError.

    This test ensures files must have BOTH valid existence AND valid extension.
    """
    class Config(luigi.Config):
        # Use an existing file that doesn't have a FASTA extension
        reference_genome = luigi.Parameter("./setup.py")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_valid():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("./tests/data/sacCer3.fa")

    assert get_genome_name(config=Config) == os.path.abspath("./tests/data/sacCer3")


def test_output_dir_passed():
    # Passed
    class Config(luigi.Config):
        output_dir = luigi.Parameter("./tests/data/")

    assert get_output_dir(config=Config) == os.path.abspath("./tests/data/")


def test_output_dir_unpassed():
    class Config(luigi.Config):
        output_dir = luigi.Parameter("")

    assert get_output_dir(config=Config) == os.getcwd()


class TestLogAndCheckCandidates:
    """Tests for log_and_check_candidates function."""

    @pytest.fixture
    def logger(self):
        """Create a test logger."""
        import logging

        return logging.getLogger("test-logger")

    def test_zero_candidates_raises_error(self, logger):
        """Zero candidates should raise ValueError with helpful hint."""
        with pytest.raises(ValueError) as exc_info:
            log_and_check_candidates(logger, "BasicFiltering", count=0)
        assert "No probes remaining" in str(exc_info.value)
        assert "BasicFiltering" in str(exc_info.value) or "TM/GC" in str(exc_info.value)

    def test_zero_candidates_stage_specific_hint(self, logger):
        """Each stage should provide a specific hint."""
        stages_and_hints = [
            ("BasicFiltering", "min-tm"),
            ("AlignProbeCandidates", "max-off-targets"),
            ("KMerFiltering", "max-kmers"),
            ("SecondaryStructureFiltering", "max-deltag"),
            ("OptimizeProbeCoverage", "spacing"),
        ]
        for stage, expected_hint in stages_and_hints:
            with pytest.raises(ValueError) as exc_info:
                log_and_check_candidates(logger, stage, count=0)
            assert expected_hint in str(exc_info.value).lower()

    def test_low_candidates_logs_warning(self, logger, capsys):
        """Less than 10 candidates should print a warning."""
        log_and_check_candidates(logger, "BasicFiltering", count=5)
        captured = capsys.readouterr().out
        assert "5" in captured and ("candidates" in captured.lower() or "remain" in captured.lower())

    def test_sufficient_candidates_no_error(self, logger):
        """10+ candidates should not raise or warn."""
        # Should not raise
        log_and_check_candidates(logger, "BasicFiltering", count=10)
        log_and_check_candidates(logger, "BasicFiltering", count=100)

    def test_count_with_previous(self, logger):
        """Previous count should be included in live display state."""
        from eFISHent.console import _steps
        _steps.clear()
        # Add a running step so print_candidate_count has something to update
        _steps.append({"name": "Test", "status": "running", "result": "", "elapsed": "", "start_time": 0})
        log_and_check_candidates(logger, "BasicFiltering", count=50, count_prev=100)
        assert "50" in _steps[0]["result"]

    def test_zero_candidates_silent_mode(self, logger):
        """Zero candidates in silent mode should still raise ValueError."""
        from unittest.mock import patch
        with patch("eFISHent.console.is_silent", return_value=True), \
             patch("eFISHent.console.print_candidate_count"), \
             patch("eFISHent.console.print_error_panel"), \
             patch("eFISHent.console.print_warning"), \
             patch("eFISHent.console.record_funnel_stage"):
            with pytest.raises(ValueError, match="No probes remaining"):
                log_and_check_candidates(logger, "BasicFiltering", count=0)

    def test_silent_mode_logs_candidates(self, logger):
        """In silent mode, candidate count should be logged, not printed."""
        from unittest.mock import patch
        with patch("eFISHent.console.is_silent", return_value=True), \
             patch("eFISHent.console.print_candidate_count") as mock_print, \
             patch("eFISHent.console.record_funnel_stage"):
            with patch.object(logger, "info") as mock_info:
                log_and_check_candidates(logger, "BasicFiltering", count=50, count_prev=100)
                mock_print.assert_not_called()
                mock_info.assert_called_once()
                call_arg = mock_info.call_args[0][0]
                assert "50" in call_arg
                assert "100" in call_arg

    def test_low_candidates_silent_mode_warns(self, logger):
        """In silent mode, low candidates should log a warning."""
        from unittest.mock import patch
        with patch("eFISHent.console.is_silent", return_value=True), \
             patch("eFISHent.console.print_candidate_count"), \
             patch("eFISHent.console.print_warning") as mock_pw, \
             patch("eFISHent.console.record_funnel_stage"):
            with patch.object(logger, "warning") as mock_warn:
                log_and_check_candidates(logger, "BasicFiltering", count=5)
                mock_pw.assert_not_called()
                mock_warn.assert_called_once()
                assert "5" in mock_warn.call_args[0][0]


class TestGetGenomeNameEdgeCases:
    """Additional edge case tests for get_genome_name."""

    def test_genome_with_fasta_extension(self):
        """Valid .fasta extension should work."""
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            f.write(b">chr1\nATCG\n")
            fasta_path = f.name

        class Config(luigi.Config):
            reference_genome = luigi.Parameter(fasta_path)

        try:
            result = get_genome_name(config=Config)
            assert result == os.path.abspath(os.path.splitext(fasta_path)[0])
        finally:
            os.unlink(fasta_path)

    def test_genome_none_raises(self):
        """None reference_genome should raise ValueError."""
        class Config(luigi.Config):
            reference_genome = luigi.Parameter(None)

        with pytest.raises(ValueError, match="Reference genome must be passed"):
            get_genome_name(config=Config)


class TestGetGeneNameEdgeCases:
    """Additional edge case tests for get_gene_name."""

    def test_no_inputs_raises(self):
        """No sequence_file, ensembl_id, or gene_name should raise ValueError."""
        class Config(luigi.Config):
            sequence_file = luigi.Parameter("")
            ensembl_id = luigi.Parameter("")
            gene_name = luigi.Parameter("")
            organism_name = luigi.Parameter("")

        with pytest.raises(ValueError, match="Could not determine gene name"):
            get_gene_name(config=Config)

    def test_gene_name_without_organism_raises(self):
        """gene_name without organism_name should raise ValueError."""
        class Config(luigi.Config):
            sequence_file = luigi.Parameter("")
            ensembl_id = luigi.Parameter("")
            gene_name = luigi.Parameter("TP53")
            organism_name = luigi.Parameter("")

        with pytest.raises(ValueError, match="Could not determine gene name"):
            get_gene_name(config=Config)


class TestGetStageDescription:
    """Tests for get_stage_description."""

    def test_known_stage(self):
        from eFISHent.util import get_stage_description
        desc = get_stage_description("PrepareSequence")
        assert desc == "Preparing gene sequence"

    def test_unknown_stage_returns_name(self):
        from eFISHent.util import get_stage_description
        desc = get_stage_description("UnknownStage")
        assert desc == "UnknownStage"


class TestHashFn:
    """Tests for hash_fn utility."""

    def test_deterministic(self):
        from eFISHent.util import hash_fn
        assert hash_fn("test") == hash_fn("test")

    def test_different_inputs(self):
        from eFISHent.util import hash_fn
        assert hash_fn("test1") != hash_fn("test2")

    def test_length(self):
        from eFISHent.util import hash_fn
        assert len(hash_fn("anything")) == 10


class TestLogStageStart:
    """Tests for log_stage_start function."""

    def test_known_stage_with_progress(self):
        """A known stage with order>0 should call print_stage."""
        import logging
        from unittest.mock import patch
        from eFISHent.util import log_stage_start

        logger = logging.getLogger("test-log-stage")
        with patch("eFISHent.config.RunConfig") as mock_rc, \
             patch("eFISHent.console.print_stage") as mock_print, \
             patch("eFISHent.console.is_silent", return_value=False):
            mock_rc.return_value.analyze_probeset = ""
            log_stage_start(logger, "PrepareSequence")
            mock_print.assert_called_once_with(1, 9, "Preparing gene sequence")

    def test_known_stage_silent_logs(self):
        """A known stage in silent mode should log info."""
        import logging
        from unittest.mock import patch

        from eFISHent.util import log_stage_start

        logger = logging.getLogger("test-log-stage-silent")
        with patch("eFISHent.config.RunConfig") as mock_rc, \
             patch("eFISHent.console.print_stage"), \
             patch("eFISHent.console.is_silent", return_value=True):
            mock_rc.return_value.analyze_probeset = ""
            with patch.object(logger, "info") as mock_info:
                log_stage_start(logger, "PrepareSequence")
                mock_info.assert_called_once_with("[1/9] Preparing gene sequence...")

    def test_index_stage_no_order(self):
        """A stage with order=0 should show as unnumbered step."""
        import logging
        from unittest.mock import patch

        from eFISHent.util import log_stage_start

        logger = logging.getLogger("test-log-stage-idx")
        with patch("eFISHent.config.RunConfig") as mock_rc, \
             patch("eFISHent.console.print_stage") as mock_print, \
             patch("eFISHent.console.is_silent", return_value=False):
            mock_rc.return_value.analyze_probeset = ""
            log_stage_start(logger, "BuildBowtieIndex")
            mock_print.assert_called_once_with(0, 0, "Building bowtie index")

    def test_unknown_stage_uses_name(self):
        """An unknown stage name should show as unnumbered step."""
        import logging
        from unittest.mock import patch

        from eFISHent.util import log_stage_start

        logger = logging.getLogger("test-log-stage-unknown")
        with patch("eFISHent.config.RunConfig") as mock_rc, \
             patch("eFISHent.console.print_stage") as mock_print, \
             patch("eFISHent.console.is_silent", return_value=False):
            mock_rc.return_value.analyze_probeset = ""
            log_stage_start(logger, "SomeNewStage")
            mock_print.assert_called_once_with(0, 0, "SomeNewStage")

    def test_analyze_mode_returns_early(self):
        """In analyze mode, log_stage_start should return immediately."""
        import logging
        from unittest.mock import patch

        from eFISHent.util import log_stage_start

        logger = logging.getLogger("test-log-stage-analyze")
        with patch("eFISHent.config.RunConfig") as mock_rc, \
             patch("eFISHent.console.print_stage") as mock_print, \
             patch("eFISHent.console.is_silent", return_value=False):
            mock_rc.return_value.analyze_probeset = "/some/probe.fasta"
            log_stage_start(logger, "PrepareSequence")
            mock_print.assert_not_called()


class TestSecureFilename:
    """Tests for secure_filename function."""

    def test_basic_filename(self):
        """Basic alphanumeric filename should pass through."""
        assert secure_filename("myfile") == "myfile"
        assert secure_filename("my_file") == "my_file"
        assert secure_filename("my-file") == "my-file"
        assert secure_filename("my.file") == "my.file"

    def test_spaces_converted(self):
        """Spaces should be converted to underscores."""
        assert secure_filename("my file") == "my_file"
        assert secure_filename("my   file") == "my_file"

    def test_special_characters_removed(self):
        """Special characters should be removed."""
        assert secure_filename("file@name") == "filename"
        assert secure_filename("file#name") == "filename"
        assert secure_filename("file$name") == "filename"
        assert secure_filename("file%name") == "filename"
        assert secure_filename("file&name") == "filename"

    def test_path_separators_removed(self):
        """Path separators should be converted to spaces then underscores."""
        result = secure_filename("path/to/file")
        assert "/" not in result
        assert "\\" not in result

    def test_unicode_normalized(self):
        """Unicode characters should be normalized to ASCII."""
        assert secure_filename("café") == "cafe"
        assert secure_filename("naïve") == "naive"
        assert secure_filename("Ångström") == "Angstrom"

    def test_leading_trailing_dots_stripped(self):
        """Leading/trailing dots and underscores should be stripped."""
        assert secure_filename(".hidden") == "hidden"
        assert secure_filename("file.") == "file"
        assert secure_filename("_private") == "private"
        assert secure_filename("file_") == "file"

    def test_empty_after_sanitization(self):
        """Filenames that become empty after sanitization."""
        assert secure_filename("...") == ""
        assert secure_filename("@#$%") == ""

    def test_real_world_examples(self):
        """Real-world filename examples."""
        assert secure_filename("Homo sapiens") == "Homo_sapiens"
        assert secure_filename("ENSG00000026025") == "ENSG00000026025"
        assert secure_filename("gene (copy)") == "gene_copy"
        assert secure_filename("résumé.pdf") == "resume.pdf"
