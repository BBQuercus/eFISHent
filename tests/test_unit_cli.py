"""Tests for CLI input validation."""

import argparse
import os
import tempfile

import pytest

from eFISHent.cli import (
    existing_count_table,
    existing_fasta_file,
    existing_gtf_file,
    non_negative_int,
    percentage,
    positive_int,
    string_to_bool,
    validate_args,
)


class TestStringToBool:
    """Tests for string_to_bool validator."""

    @pytest.mark.parametrize(
        "value", ["true", "True", "TRUE", "yes", "Yes", "y", "1", "t"]
    )
    def test_truthy_values(self, value):
        assert string_to_bool(value) is True

    @pytest.mark.parametrize(
        "value", ["false", "False", "FALSE", "no", "No", "n", "0", "f"]
    )
    def test_falsy_values(self, value):
        assert string_to_bool(value) is False

    def test_bool_passthrough(self):
        assert string_to_bool(True) is True
        assert string_to_bool(False) is False

    @pytest.mark.parametrize("value", ["maybe", "2", "yep", "nope", ""])
    def test_invalid_values(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            string_to_bool(value)
        assert "Invalid boolean value" in str(exc_info.value)


class TestPositiveInt:
    """Tests for positive_int validator."""

    @pytest.mark.parametrize("value,expected", [("1", 1), ("5", 5), ("100", 100)])
    def test_valid_values(self, value, expected):
        assert positive_int(value) == expected

    @pytest.mark.parametrize("value", ["0", "-1", "-100"])
    def test_zero_and_negative(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            positive_int(value)
        assert "must be >= 1" in str(exc_info.value)

    @pytest.mark.parametrize("value", ["abc", "1.5", "", "one"])
    def test_non_integer(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            positive_int(value)
        assert "Invalid integer" in str(exc_info.value)


class TestNonNegativeInt:
    """Tests for non_negative_int validator."""

    @pytest.mark.parametrize("value,expected", [("0", 0), ("1", 1), ("100", 100)])
    def test_valid_values(self, value, expected):
        assert non_negative_int(value) == expected

    @pytest.mark.parametrize("value", ["-1", "-100"])
    def test_negative(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            non_negative_int(value)
        assert "must be >= 0" in str(exc_info.value)

    @pytest.mark.parametrize("value", ["abc", "1.5", ""])
    def test_non_integer(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            non_negative_int(value)
        assert "Invalid integer" in str(exc_info.value)


class TestPercentage:
    """Tests for percentage validator."""

    @pytest.mark.parametrize(
        "value,expected", [("0", 0.0), ("50", 50.0), ("100", 100.0), ("33.5", 33.5)]
    )
    def test_valid_values(self, value, expected):
        assert percentage(value) == expected

    @pytest.mark.parametrize("value", ["-1", "-0.1", "100.1", "150", "200"])
    def test_out_of_range(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            percentage(value)
        assert "must be 0-100" in str(exc_info.value)

    @pytest.mark.parametrize("value", ["abc", "", "fifty"])
    def test_non_numeric(self, value):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            percentage(value)
        assert "Invalid number" in str(exc_info.value)


class TestExistingFastaFile:
    """Tests for existing_fasta_file validator."""

    def test_empty_value_allowed(self):
        assert existing_fasta_file("") == ""

    def test_valid_fasta_file(self):
        result = existing_fasta_file("./tests/data/sacCer3.fa")
        assert result == "./tests/data/sacCer3.fa"

    def test_nonexistent_file(self):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            existing_fasta_file("/nonexistent/path/genome.fa")
        assert "File not found" in str(exc_info.value)

    def test_wrong_extension(self):
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"test")
            temp_path = f.name
        try:
            with pytest.raises(argparse.ArgumentTypeError) as exc_info:
                existing_fasta_file(temp_path)
            assert "Expected FASTA file" in str(exc_info.value)
        finally:
            os.unlink(temp_path)


class TestExistingGtfFile:
    """Tests for existing_gtf_file validator."""

    def test_empty_value_allowed(self):
        assert existing_gtf_file("") == ""

    def test_nonexistent_file(self):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            existing_gtf_file("/nonexistent/path/annotation.gtf")
        assert "File not found" in str(exc_info.value)

    def test_wrong_extension(self):
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"test")
            temp_path = f.name
        try:
            with pytest.raises(argparse.ArgumentTypeError) as exc_info:
                existing_gtf_file(temp_path)
            assert "Expected GTF/GFF file" in str(exc_info.value)
        finally:
            os.unlink(temp_path)


class TestExistingCountTable:
    """Tests for existing_count_table validator."""

    def test_empty_value_allowed(self):
        assert existing_count_table("") == ""

    def test_nonexistent_file(self):
        with pytest.raises(argparse.ArgumentTypeError) as exc_info:
            existing_count_table("/nonexistent/path/counts.csv")
        assert "File not found" in str(exc_info.value)

    def test_wrong_extension(self):
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            f.write(b"test")
            temp_path = f.name
        try:
            with pytest.raises(argparse.ArgumentTypeError) as exc_info:
                existing_count_table(temp_path)
            assert "Expected count table" in str(exc_info.value)
        finally:
            os.unlink(temp_path)

    @pytest.mark.parametrize("suffix", [".csv", ".tsv", ".txt"])
    def test_valid_extensions(self, suffix):
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
            f.write(b"gene_id,count\nENSG001,100")
            temp_path = f.name
        try:
            result = existing_count_table(temp_path)
            assert result == temp_path
        finally:
            os.unlink(temp_path)


class TestValidateArgs:
    """Tests for validate_args cross-parameter validation."""

    @pytest.fixture
    def mock_parser(self):
        """Create a mock parser that captures errors."""
        parser = argparse.ArgumentParser()
        return parser

    @pytest.fixture
    def default_args(self):
        """Create default valid args namespace."""
        return argparse.Namespace(
            # General
            reference_genome="./tests/data/sacCer3.fa",
            reference_annotation="",
            threads=2,
            output_dir="",
            # Run
            build_indices=False,
            analyze_probeset="",
            save_intermediates=False,
            optimization_method="greedy",
            optimization_time_limit=60,
            # Sequence
            sequence_file="./tests/data/renilla.fasta",
            ensembl_id="",
            gene_name="",
            organism_name="",
            is_plus_strand=True,
            is_endogenous=True,
            # Probe
            min_length=21,
            max_length=25,
            spacing=2,
            min_tm=40.0,
            max_tm=60.0,
            min_gc=20.0,
            max_gc=80.0,
            formamide_concentration=10.0,
            na_concentration=330.0,
            max_off_targets=0,
            no_alternative_loci=False,
            encode_count_table="",
            max_expression_percentage=50.0,
            kmer_length=15,
            max_kmers=5,
            max_deltag=-10.0,
            sequence_similarity=0,
        )

    def test_valid_args_pass(self, mock_parser, default_args):
        """Valid args should not raise."""
        validate_args(default_args, mock_parser)  # Should not raise

    def test_min_length_greater_than_max_length(self, mock_parser, default_args):
        default_args.min_length = 30
        default_args.max_length = 20
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_min_tm_greater_than_max_tm(self, mock_parser, default_args):
        default_args.min_tm = 70.0
        default_args.max_tm = 50.0
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_min_gc_greater_than_max_gc(self, mock_parser, default_args):
        default_args.min_gc = 90.0
        default_args.max_gc = 10.0
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_kmer_length_too_large(self, mock_parser, default_args):
        default_args.kmer_length = 25  # >= min_length of 21
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_no_sequence_input(self, mock_parser, default_args):
        default_args.sequence_file = ""
        default_args.ensembl_id = ""
        default_args.gene_name = ""
        default_args.organism_name = ""
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_gene_name_without_organism(self, mock_parser, default_args):
        default_args.sequence_file = ""
        default_args.gene_name = "ACTB"
        default_args.organism_name = ""
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_organism_without_gene_name(self, mock_parser, default_args):
        default_args.sequence_file = ""
        default_args.gene_name = ""
        default_args.organism_name = "Homo sapiens"
        with pytest.raises(SystemExit):
            validate_args(default_args, mock_parser)

    def test_gene_name_with_organism_valid(self, mock_parser, default_args):
        default_args.sequence_file = ""
        default_args.gene_name = "ACTB"
        default_args.organism_name = "Homo sapiens"
        validate_args(default_args, mock_parser)  # Should not raise

    def test_ensembl_id_valid(self, mock_parser, default_args):
        default_args.sequence_file = ""
        default_args.ensembl_id = "ENSG00000001234"
        validate_args(default_args, mock_parser)  # Should not raise

    def test_encode_count_table_without_annotation(self, mock_parser, default_args):
        # Create a temp count table file
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            f.write(b"gene_id,count\nENSG001,100")
            temp_path = f.name
        try:
            default_args.encode_count_table = temp_path
            default_args.reference_annotation = ""
            with pytest.raises(SystemExit):
                validate_args(default_args, mock_parser)
        finally:
            os.unlink(temp_path)

    def test_build_indices_skips_sequence_check(self, mock_parser, default_args):
        default_args.build_indices = True
        default_args.sequence_file = ""
        default_args.ensembl_id = ""
        default_args.gene_name = ""
        default_args.organism_name = ""
        validate_args(default_args, mock_parser)  # Should not raise

    def test_analyze_probeset_skips_sequence_check(self, mock_parser, default_args):
        default_args.analyze_probeset = "./tests/data/renilla.fasta"
        default_args.sequence_file = ""
        default_args.ensembl_id = ""
        default_args.gene_name = ""
        default_args.organism_name = ""
        validate_args(default_args, mock_parser)  # Should not raise
