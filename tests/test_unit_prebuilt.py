"""Tests for pre-built genome index management."""

import json
import os
import tempfile

import pytest

from eFISHent.prebuilt import (
    GENOME_ALIASES,
    get_cache_dir,
    get_genome_dir,
    get_reference_paths,
    is_genome_cached,
    list_available_genomes,
    resolve_genome,
    verify_checksums,
)


class TestResolveGenome:
    def test_resolve_known_aliases(self):
        assert resolve_genome("hg38") == "homo_sapiens/GRCh38"
        assert resolve_genome("human") == "homo_sapiens/GRCh38"
        assert resolve_genome("GRCh38") == "homo_sapiens/GRCh38"
        assert resolve_genome("mm39") == "mus_musculus/GRCm39"
        assert resolve_genome("mouse") == "mus_musculus/GRCm39"
        assert resolve_genome("dm6") == "drosophila_melanogaster/BDGP6"
        assert resolve_genome("fly") == "drosophila_melanogaster/BDGP6"

    def test_resolve_unknown_alias(self):
        with pytest.raises(ValueError, match="Unknown genome"):
            resolve_genome("invalid_genome")

    def test_all_aliases_resolve(self):
        """Every alias in GENOME_ALIASES should resolve without error."""
        for alias in GENOME_ALIASES:
            result = resolve_genome(alias)
            assert isinstance(result, str)
            assert "/" in result


class TestCacheManagement:
    def test_get_cache_dir_creates_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = os.path.join(tmpdir, "test_cache")
            result = get_cache_dir(cache)
            assert os.path.isdir(result)
            assert result == cache

    def test_get_genome_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = get_genome_dir("homo_sapiens/GRCh38", tmpdir)
            assert result == os.path.join(tmpdir, "homo_sapiens/GRCh38")

    def test_is_genome_cached_false_when_empty(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            assert is_genome_cached("homo_sapiens/GRCh38", tmpdir) is False

    def test_is_genome_cached_true_with_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            # Create metadata with one file
            dummy_file = os.path.join(genome_dir, "genome.fa")
            with open(dummy_file, "w") as f:
                f.write(">chr1\nATCG\n")
            metadata = {"files": {"genome.fa": "abc123"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            assert is_genome_cached("homo_sapiens/GRCh38", tmpdir) is True

    def test_is_genome_cached_false_missing_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            metadata = {"files": {"genome.fa": "abc123"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            # genome.fa not created
            assert is_genome_cached("homo_sapiens/GRCh38", tmpdir) is False


class TestGetReferencePaths:
    def test_raises_when_not_cached(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(FileNotFoundError, match="not cached"):
                get_reference_paths("homo_sapiens/GRCh38", tmpdir)

    def test_returns_paths_when_cached(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            for name in ["genome.fa", "annotation.gtf", "transcriptome.fa"]:
                with open(os.path.join(genome_dir, name), "w") as f:
                    f.write("dummy")
            metadata = {
                "files": {
                    "genome.fa": "x",
                    "annotation.gtf": "x",
                    "transcriptome.fa": "x",
                }
            }
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)

            paths = get_reference_paths("homo_sapiens/GRCh38", tmpdir)
            assert "reference_genome" in paths
            assert paths["reference_genome"].endswith("genome.fa")
            assert "reference_annotation" in paths
            assert "reference_transcriptome" in paths


class TestListAvailableGenomes:
    def test_returns_list(self):
        genomes = list_available_genomes()
        assert len(genomes) >= 3  # human, mouse, fly

    def test_each_entry_has_required_fields(self):
        genomes = list_available_genomes()
        for genome in genomes:
            assert "id" in genome
            assert "aliases" in genome
            assert "display" in genome
            assert "/" in genome["id"]


class TestVerifyChecksums:
    def test_missing_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            valid, failed = verify_checksums("homo_sapiens/GRCh38", tmpdir)
            assert valid is False
            assert "metadata.json missing" in failed

    def test_valid_checksum(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "test/genome")
            os.makedirs(genome_dir)
            # Create a file and compute its hash
            content = b"test content"
            filepath = os.path.join(genome_dir, "test.fa")
            with open(filepath, "wb") as f:
                f.write(content)
            expected_hash = hashlib.sha256(content).hexdigest()
            metadata = {"files": {"test.fa": expected_hash}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)

            valid, failed = verify_checksums("test/genome", tmpdir)
            assert valid is True
            assert failed == []
