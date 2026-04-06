"""Tests for pre-built genome index management."""

import json
import os
import tempfile

import pytest

from eFISHent.prebuilt import (
    GENOME_ALIASES,
    PLANNED_GENOMES,
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

    def test_resolve_new_aliases(self):
        """Yeast and C. elegans should resolve (indices uploaded)."""
        assert resolve_genome("yeast") == "saccharomyces_cerevisiae/R64"
        assert resolve_genome("sacCer3") == "saccharomyces_cerevisiae/R64"
        assert resolve_genome("R64") == "saccharomyces_cerevisiae/R64"
        assert resolve_genome("worm") == "caenorhabditis_elegans/WBcel235"
        assert resolve_genome("ce11") == "caenorhabditis_elegans/WBcel235"
        assert resolve_genome("elegans") == "caenorhabditis_elegans/WBcel235"

    def test_resolve_planned_genome_raises(self):
        """Planned genomes should raise with a helpful 'not yet available' message."""
        with pytest.raises(ValueError, match="not yet available"):
            resolve_genome("dm6")
        with pytest.raises(ValueError, match="not yet available"):
            resolve_genome("fly")
        with pytest.raises(ValueError, match="Zebrafish"):
            resolve_genome("danRer11")
        with pytest.raises(ValueError, match="Rat"):
            resolve_genome("rn7")

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
        assert len(genomes) >= 2  # human, mouse

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

    def test_corrupt_metadata_json(self):
        """Corrupt metadata.json should return False with 'corrupt' message."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "test/genome")
            os.makedirs(genome_dir)
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                f.write("{invalid json")
            valid, failed = verify_checksums("test/genome", tmpdir)
            assert valid is False
            assert "metadata.json corrupt" in failed

    def test_missing_file_in_checksums(self):
        """A file listed in metadata but not on disk should appear in failed."""
        import hashlib
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "test/genome")
            os.makedirs(genome_dir)
            metadata = {"files": {"missing.fa": "abc123", "also_missing.bt2": "def456"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            valid, failed = verify_checksums("test/genome", tmpdir)
            assert valid is False
            assert "missing.fa" in failed
            assert "also_missing.bt2" in failed

    def test_invalid_checksum(self):
        """A file with wrong checksum should appear in failed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "test/genome")
            os.makedirs(genome_dir)
            filepath = os.path.join(genome_dir, "test.fa")
            with open(filepath, "wb") as f:
                f.write(b"test content")
            metadata = {"files": {"test.fa": "0000000000000000000000000000000000000000000000000000000000000000"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            valid, failed = verify_checksums("test/genome", tmpdir)
            assert valid is False
            assert "test.fa" in failed

    def test_empty_files_dict_is_valid(self):
        """metadata with empty files dict should be considered valid."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "test/genome")
            os.makedirs(genome_dir)
            metadata = {"files": {}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            valid, failed = verify_checksums("test/genome", tmpdir)
            assert valid is True
            assert failed == []


class TestIsGenomeCachedEdgeCases:
    def test_corrupt_metadata_returns_false(self):
        """Corrupt JSON in metadata.json should return False."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                f.write("not valid json{{{")
            assert is_genome_cached("homo_sapiens/GRCh38", tmpdir) is False

    def test_metadata_without_files_key(self):
        """metadata.json without 'files' key should still return True (no files to check)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            metadata = {"version": "1.0"}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)
            # .get("files", {}) returns {}, so no files to check -> True
            assert is_genome_cached("homo_sapiens/GRCh38", tmpdir) is True


class TestDownloadGenome:
    def test_download_already_cached_skips(self):
        """download_genome should return early when already cached."""
        from unittest.mock import patch
        from eFISHent.prebuilt import download_genome

        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            dummy_file = os.path.join(genome_dir, "genome.fa")
            with open(dummy_file, "w") as f:
                f.write(">chr1\nATCG\n")
            metadata = {"files": {"genome.fa": "abc123"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)

            # Should not attempt to import huggingface_hub
            result = download_genome("homo_sapiens/GRCh38", cache_dir=tmpdir)
            assert result == genome_dir

    def test_download_force_ignores_cache(self):
        """download_genome with force=True should not skip even if cached."""
        from unittest.mock import patch, MagicMock
        from eFISHent.prebuilt import download_genome

        with tempfile.TemporaryDirectory() as tmpdir:
            genome_dir = os.path.join(tmpdir, "homo_sapiens/GRCh38")
            os.makedirs(genome_dir)
            dummy_file = os.path.join(genome_dir, "genome.fa")
            with open(dummy_file, "w") as f:
                f.write(">chr1\nATCG\n")
            metadata = {"files": {"genome.fa": "abc123"}}
            with open(os.path.join(genome_dir, "metadata.json"), "w") as f:
                json.dump(metadata, f)

            mock_item = MagicMock()
            mock_item.rfilename = "homo_sapiens/GRCh38/genome.fa"

            with patch.dict("sys.modules", {"huggingface_hub": MagicMock()}):
                import sys
                hf_mock = sys.modules["huggingface_hub"]
                hf_mock.list_repo_tree.return_value = [mock_item]
                hf_mock.hf_hub_download.return_value = dummy_file

                with patch("eFISHent.console.console", MagicMock()):
                    result = download_genome(
                        "homo_sapiens/GRCh38", cache_dir=tmpdir, force=True
                    )
            assert result == genome_dir

    def test_download_no_files_raises(self):
        """download_genome should raise when HF returns no files."""
        from unittest.mock import patch, MagicMock
        from eFISHent.prebuilt import download_genome

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict("sys.modules", {"huggingface_hub": MagicMock()}):
                import sys
                hf_mock = sys.modules["huggingface_hub"]
                hf_mock.list_repo_tree.return_value = []

                with patch("eFISHent.console.console", MagicMock()):
                    with pytest.raises(RuntimeError, match="No files found"):
                        download_genome("homo_sapiens/GRCh38", cache_dir=tmpdir)

    def test_download_list_error_raises(self):
        """download_genome should raise when listing files fails."""
        from unittest.mock import patch, MagicMock
        from eFISHent.prebuilt import download_genome

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict("sys.modules", {"huggingface_hub": MagicMock()}):
                import sys
                hf_mock = sys.modules["huggingface_hub"]
                hf_mock.list_repo_tree.side_effect = Exception("Network error")

                with patch("eFISHent.console.console", MagicMock()):
                    with pytest.raises(RuntimeError, match="Could not list files"):
                        download_genome("homo_sapiens/GRCh38", cache_dir=tmpdir)


class TestResolveGenomeErrorMessage:
    def test_error_lists_available_genomes(self):
        """Error message for truly unknown alias should list available genomes."""
        with pytest.raises(ValueError) as exc_info:
            resolve_genome("notareal_genome")
        msg = str(exc_info.value)
        assert "homo_sapiens/GRCh38" in msg
        assert "mus_musculus/GRCm39" in msg
        assert "Aliases:" in msg

    def test_planned_genome_error_is_helpful(self):
        """Planned genome error should name the organism and link to wiki."""
        with pytest.raises(ValueError) as exc_info:
            resolve_genome("zebrafish")
        msg = str(exc_info.value)
        assert "Zebrafish" in msg
        assert "not yet available" in msg
        assert "wiki" in msg.lower() or "instructions" in msg.lower()

    def test_error_lists_aliases(self):
        """Error message should list available aliases."""
        with pytest.raises(ValueError) as exc_info:
            resolve_genome("notreal")
        msg = str(exc_info.value)
        assert "hg38" in msg
        assert "mm39" in msg


class TestListAvailableGenomesDetails:
    def test_aliases_are_sorted(self):
        genomes = list_available_genomes()
        for genome in genomes:
            assert genome["aliases"] == sorted(genome["aliases"])

    def test_display_format(self):
        genomes = list_available_genomes()
        for genome in genomes:
            assert genome["id"] in genome["display"]
            assert "(" in genome["display"]
            assert ")" in genome["display"]

    def test_all_canonical_ids_covered(self):
        """All unique canonical IDs should appear."""
        genomes = list_available_genomes()
        canonical_ids = {g["id"] for g in genomes}
        expected_ids = set(GENOME_ALIASES.values())
        assert canonical_ids == expected_ids
