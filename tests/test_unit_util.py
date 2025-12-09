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

    def test_low_candidates_logs_warning(self, logger, caplog):
        """Less than 10 candidates should log a warning."""
        import logging

        with caplog.at_level(logging.WARNING):
            log_and_check_candidates(logger, "BasicFiltering", count=5)
        assert "5 candidates remain" in caplog.text or "Only 5" in caplog.text

    def test_sufficient_candidates_no_error(self, logger):
        """10+ candidates should not raise or warn."""
        # Should not raise
        log_and_check_candidates(logger, "BasicFiltering", count=10)
        log_and_check_candidates(logger, "BasicFiltering", count=100)

    def test_count_with_previous(self, logger, caplog):
        """Previous count should be included in log message."""
        import logging

        with caplog.at_level(logging.INFO):
            log_and_check_candidates(logger, "BasicFiltering", count=50, count_prev=100)
        assert "50" in caplog.text
        assert "100" in caplog.text


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
