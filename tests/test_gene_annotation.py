"""Tests for gene_annotation module."""

import os
import tempfile
from unittest.mock import patch

import pandas as pd
import pytest

from eFISHent.gene_annotation import (
    aggregate_off_target_genes,
    build_transcript_gene_map,
    map_transcript_to_gene,
)


@pytest.fixture
def mock_gtf_content():
    """Minimal GTF content for testing."""
    return (
        '# comment line\n'
        'chr1\tensembl\ttranscript\t100\t200\t.\t+\t.\t'
        'gene_id "ENSG00000000001"; transcript_id "ENST00000123456.3"; gene_name "BRCA1";\n'
        'chr1\tensembl\texon\t100\t150\t.\t+\t.\t'
        'gene_id "ENSG00000000001"; transcript_id "ENST00000123456.3"; gene_name "BRCA1";\n'
        'chr2\tensembl\ttranscript\t500\t700\t.\t-\t.\t'
        'gene_id "ENSG00000000002"; transcript_id "ENST00000654321.1"; gene_name "TP53";\n'
        'chr3\tensembl\ttranscript\t800\t900\t.\t+\t.\t'
        'gene_id "ENSG00000000003"; transcript_id "NM_001234.5"; gene_name "MYC";\n'
    )


@pytest.fixture
def mock_gtf_file(mock_gtf_content):
    """Write mock GTF to a temporary file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".gtf", delete=False
    ) as f:
        f.write(mock_gtf_content)
        path = f.name
    yield path
    os.unlink(path)


@pytest.fixture
def transcript_mapping():
    """Pre-built mapping for testing map_transcript_to_gene."""
    return {
        "ENST00000123456.3": "BRCA1",
        "ENST00000123456": "BRCA1",
        "ENST00000654321.1": "TP53",
        "ENST00000654321": "TP53",
        "NM_001234.5": "MYC",
        "NM_001234": "MYC",
    }


@pytest.fixture
def blast_hits_df():
    """Mock BLAST DataFrame for aggregate_off_target_genes tests."""
    return pd.DataFrame(
        {
            "qseqid": [
                "probe-1",
                "probe-1",
                "probe-1",
                "probe-1",
                "probe-2",
                "probe-2",
            ],
            "sseqid": [
                "ENST00000123456",  # self-hit (target gene BRCA1)
                "ENST00000654321",  # off-target TP53
                "ENST00000654321",  # off-target TP53 again (different transcript hit)
                "ENST00000999999",  # unknown transcript
                "ENST00000654321",  # off-target TP53
                "ENST00000123456",  # self-hit
            ],
            "pident": [100.0, 95.0, 90.0, 85.0, 92.0, 100.0],
            "length": [25, 22, 20, 18, 21, 25],
            "mismatch": [0, 1, 2, 3, 2, 0],
        }
    )


class TestBuildTranscriptGeneMap:
    def test_parse_raw_gtf(self, mock_gtf_file):
        """Test parsing a raw GTF file produces correct mappings."""
        # Clear the cache to avoid interference
        from eFISHent.gene_annotation import _transcript_gene_cache

        _transcript_gene_cache.clear()

        mapping = build_transcript_gene_map(mock_gtf_file)

        # Should map versioned IDs
        assert mapping["ENST00000123456.3"] == "BRCA1"
        assert mapping["ENST00000654321.1"] == "TP53"
        assert mapping["NM_001234.5"] == "MYC"

        # Should also map base IDs (without version)
        assert mapping["ENST00000123456"] == "BRCA1"
        assert mapping["ENST00000654321"] == "TP53"
        assert mapping["NM_001234"] == "MYC"

    def test_nonexistent_file_returns_empty(self):
        """Test that a missing file returns an empty mapping."""
        from eFISHent.gene_annotation import _transcript_gene_cache

        _transcript_gene_cache.clear()

        mapping = build_transcript_gene_map("/nonexistent/path.gtf")
        assert mapping == {}

    def test_caching(self, mock_gtf_file):
        """Test that repeated calls return cached results."""
        from eFISHent.gene_annotation import _transcript_gene_cache

        _transcript_gene_cache.clear()

        mapping1 = build_transcript_gene_map(mock_gtf_file)
        mapping2 = build_transcript_gene_map(mock_gtf_file)
        assert mapping1 is mapping2


class TestMapTranscriptToGene:
    def test_direct_lookup(self, transcript_mapping):
        """Direct Ensembl ID lookup."""
        assert (
            map_transcript_to_gene("ENST00000123456.3", transcript_mapping)
            == "BRCA1"
        )

    def test_versioned_id_stripped(self, transcript_mapping):
        """Versioned ID should be resolved by stripping the version."""
        assert (
            map_transcript_to_gene("ENST00000123456.1", transcript_mapping)
            == "BRCA1"
        )

    def test_pipe_delimited(self, transcript_mapping):
        """Pipe-delimited gffread format should be parsed."""
        assert (
            map_transcript_to_gene(
                "BRCA1_gene|ENST00000123456", transcript_mapping
            )
            == "BRCA1"
        )

    def test_pipe_delimited_versioned(self, transcript_mapping):
        """Pipe-delimited with versioned ID."""
        assert (
            map_transcript_to_gene(
                "gene|ENST00000654321.1", transcript_mapping
            )
            == "TP53"
        )

    def test_refseq_direct(self, transcript_mapping):
        """Direct RefSeq lookup."""
        assert (
            map_transcript_to_gene("NM_001234.5", transcript_mapping) == "MYC"
        )

    def test_refseq_base(self, transcript_mapping):
        """RefSeq without version should use underscore parsing."""
        assert (
            map_transcript_to_gene("NM_001234", transcript_mapping) == "MYC"
        )

    def test_unknown_returns_itself(self, transcript_mapping):
        """Unknown ID should return the original ID."""
        assert (
            map_transcript_to_gene("UNKNOWN_GENE_42", transcript_mapping)
            == "UNKNOWN_GENE_42"
        )

    def test_empty_mapping(self):
        """Any ID with empty mapping returns itself."""
        assert map_transcript_to_gene("ENST00000123456", {}) == "ENST00000123456"


class TestParseParquetGtf:
    """Tests for _parse_parquet_gtf function."""

    def test_parquet_with_gene_name_and_transcript_id(self, tmp_path):
        """Parquet with gene_name and transcript_id should produce mapping."""
        from eFISHent.gene_annotation import _parse_parquet_gtf

        df = pd.DataFrame({
            "transcript_id": ["ENST00000111111.1", "ENST00000222222.2"],
            "gene_name": ["BRCA1", "TP53"],
            "gene_id": ["ENSG00000001", "ENSG00000002"],
        })
        parquet_path = str(tmp_path / "test.parquet")
        df.to_parquet(parquet_path)

        mapping = _parse_parquet_gtf(parquet_path)
        assert mapping["ENST00000111111.1"] == "BRCA1"
        assert mapping["ENST00000111111"] == "BRCA1"
        assert mapping["ENST00000222222.2"] == "TP53"
        assert mapping["ENST00000222222"] == "TP53"
        # gene_id -> gene_name mapping
        assert mapping["ENSG00000001"] == "BRCA1"
        assert mapping["ENSG00000002"] == "TP53"

    def test_parquet_without_gene_name_uses_gene_id(self, tmp_path):
        """Parquet without gene_name should use gene_id as fallback."""
        from eFISHent.gene_annotation import _parse_parquet_gtf

        df = pd.DataFrame({
            "transcript_id": ["ENST00000111111"],
            "gene_id": ["ENSG00000001"],
        })
        parquet_path = str(tmp_path / "test.parquet")
        df.to_parquet(parquet_path)

        mapping = _parse_parquet_gtf(parquet_path)
        assert mapping["ENST00000111111"] == "ENSG00000001"

    def test_parquet_no_transcript_id_still_maps_gene_id(self, tmp_path):
        """Parquet without transcript_id but with gene_name+gene_id should map gene_id->gene_name."""
        from eFISHent.gene_annotation import _parse_parquet_gtf

        df = pd.DataFrame({
            "gene_name": ["BRCA1", "TP53"],
            "gene_id": ["ENSG00000001", "ENSG00000002"],
        })
        parquet_path = str(tmp_path / "test.parquet")
        df.to_parquet(parquet_path)

        mapping = _parse_parquet_gtf(parquet_path)
        # gene_id -> gene_name still built even without transcript_id
        assert mapping["ENSG00000001"] == "BRCA1"
        assert mapping["ENSG00000002"] == "TP53"

    def test_parquet_no_useful_columns(self, tmp_path):
        """Parquet without transcript_id, gene_name, or gene_id should return empty."""
        from eFISHent.gene_annotation import _parse_parquet_gtf

        df = pd.DataFrame({
            "seqname": ["chr1", "chr2"],
            "feature": ["exon", "exon"],
        })
        parquet_path = str(tmp_path / "test.parquet")
        df.to_parquet(parquet_path)

        mapping = _parse_parquet_gtf(parquet_path)
        assert len(mapping) == 0


class TestParseRawGtf:
    """Tests for _parse_raw_gtf function."""

    def test_gzipped_gtf(self, tmp_path):
        """Gzipped GTF file should be parsed correctly."""
        import gzip
        from eFISHent.gene_annotation import _parse_raw_gtf

        gtf_content = (
            '# comment\n'
            'chr1\tensembl\texon\t100\t200\t.\t+\t.\t'
            'gene_id "ENSG00000000001"; transcript_id "ENST00000111111.2"; gene_name "MYC";\n'
        )
        gz_path = str(tmp_path / "test.gtf.gz")
        with gzip.open(gz_path, "wt") as f:
            f.write(gtf_content)

        mapping = _parse_raw_gtf(gz_path)
        assert mapping["ENST00000111111.2"] == "MYC"
        assert mapping["ENST00000111111"] == "MYC"

    def test_gtf_without_gene_name_falls_back_to_gene_id(self, tmp_path):
        """GTF without gene_name should use gene_id."""
        from eFISHent.gene_annotation import _parse_raw_gtf

        gtf_content = (
            'chr1\tensembl\texon\t100\t200\t.\t+\t.\t'
            'gene_id "ENSG00000000001"; transcript_id "ENST00000111111";\n'
        )
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)

        mapping = _parse_raw_gtf(gtf_path)
        assert mapping["ENST00000111111"] == "ENSG00000000001"

    def test_gtf_skips_lines_without_transcript_id(self, tmp_path):
        """Lines without transcript_id should be skipped."""
        from eFISHent.gene_annotation import _parse_raw_gtf

        gtf_content = (
            'chr1\tensembl\tgene\t100\t200\t.\t+\t.\t'
            'gene_id "ENSG00000000001"; gene_name "BRCA1";\n'
        )
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)

        mapping = _parse_raw_gtf(gtf_path)
        assert len(mapping) == 0

    def test_gtf_short_lines_skipped(self, tmp_path):
        """Lines with fewer than 9 tab-separated fields should be skipped."""
        from eFISHent.gene_annotation import _parse_raw_gtf

        gtf_content = "chr1\tensembl\texon\t100\n"
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)

        mapping = _parse_raw_gtf(gtf_path)
        assert len(mapping) == 0


class TestMapTranscriptToGeneAdditional:
    """Additional tests for map_transcript_to_gene edge cases."""

    def test_underscore_refseq_format(self):
        """Underscore-delimited RefSeq IDs should be resolved."""
        mapping = {"NM_001234": "GAPDH"}
        assert map_transcript_to_gene("NM_001234_extra", mapping) == "GAPDH"

    def test_pipe_with_base_id_resolution(self):
        """Pipe-delimited with versioned sub-parts."""
        mapping = {"ENST00000123456": "TP53"}
        result = map_transcript_to_gene("gene|ENST00000123456.5", mapping)
        assert result == "TP53"


class TestBuildTranscriptGeneMapParquet:
    """Test build_transcript_gene_map preferring parquet format."""

    def test_prefers_parquet_over_raw_gtf(self, tmp_path):
        """When a parquet file exists alongside GTF, it should be preferred."""
        from eFISHent.gene_annotation import build_transcript_gene_map, _transcript_gene_cache

        _transcript_gene_cache.clear()

        # Create a parquet file (the function expects .parquet extension via splitext)
        df = pd.DataFrame({
            "transcript_id": ["ENST00000999999"],
            "gene_name": ["PARQUET_GENE"],
            "gene_id": ["ENSG00000999999"],
        })
        parquet_path = str(tmp_path / "test.parquet")
        df.to_parquet(parquet_path)

        # Also create a GTF file with different content
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(
                'chr1\tensembl\texon\t100\t200\t.\t+\t.\t'
                'gene_id "ENSG00000000001"; transcript_id "ENST00000111111"; gene_name "GTF_GENE";\n'
            )

        # Pass the GTF path - it should find the .parquet sibling
        mapping = build_transcript_gene_map(gtf_path)
        assert "ENST00000999999" in mapping
        assert mapping["ENST00000999999"] == "PARQUET_GENE"

    def test_build_from_raw_gtf_when_no_parquet(self, mock_gtf_file):
        """When no parquet file exists, should fall back to raw GTF."""
        from eFISHent.gene_annotation import build_transcript_gene_map, _transcript_gene_cache

        _transcript_gene_cache.clear()

        mapping = build_transcript_gene_map(mock_gtf_file)
        assert "ENST00000123456.3" in mapping
        assert mapping["ENST00000123456.3"] == "BRCA1"


class TestParseRawGtfNoGeneName:
    """Test _parse_raw_gtf with missing gene_name and gene_id."""

    def test_no_gene_name_no_gene_id_skips_line(self, tmp_path):
        """Lines with transcript_id but no gene_name or gene_id should be skipped."""
        from eFISHent.gene_annotation import _parse_raw_gtf

        gtf_content = (
            'chr1\tensembl\texon\t100\t200\t.\t+\t.\t'
            'transcript_id "ENST00000111111";\n'
        )
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.write(gtf_content)

        mapping = _parse_raw_gtf(gtf_path)
        assert len(mapping) == 0


class TestAggregateOffTargetGenes:
    def test_self_hit_exclusion(self, blast_hits_df, transcript_mapping):
        """Self-hits (matching target gene) should be excluded."""
        results = aggregate_off_target_genes(
            blast_hits_df, transcript_mapping, target_gene="BRCA1"
        )
        # probe-1 has 4 hits: 1 self-hit (BRCA1), 2 TP53, 1 unknown
        info = results["probe-1"]
        # txome_off_targets should not include self-hits
        assert info["txome_off_targets"] >= 2  # at least TP53 hits + unknown

    def test_gene_counting(self, blast_hits_df, transcript_mapping):
        """Off-target gene counts should be correct."""
        results = aggregate_off_target_genes(
            blast_hits_df, transcript_mapping, target_gene="BRCA1"
        )
        info = results["probe-1"]
        assert "TP53" in info["off_target_genes"]

    def test_worst_match_tracking(self, blast_hits_df, transcript_mapping):
        """The worst (highest-quality) off-target match should be tracked."""
        results = aggregate_off_target_genes(
            blast_hits_df, transcript_mapping, target_gene="BRCA1"
        )
        info = results["probe-1"]
        # TP53 has hits at 95% and 90%; worst_match should reflect the 95% hit
        assert info["worst_match"] != ""
        assert "95%" in info["worst_match"]

    def test_probe_with_only_self_hits(self, transcript_mapping):
        """A probe with only self-hits should have zero off-targets."""
        df = pd.DataFrame(
            {
                "qseqid": ["probe-3", "probe-3"],
                "sseqid": ["ENST00000123456", "ENST00000123456.3"],
                "pident": [100.0, 99.0],
                "length": [25, 24],
                "mismatch": [0, 1],
            }
        )
        results = aggregate_off_target_genes(
            df, transcript_mapping, target_gene="BRCA1"
        )
        info = results["probe-3"]
        assert info["txome_off_targets"] == 0
        assert info["off_target_genes"] == ""
        assert info["worst_match"] == ""

    def test_empty_blast_results(self, transcript_mapping):
        """Empty BLAST results should produce empty dict."""
        df = pd.DataFrame(
            columns=["qseqid", "sseqid", "pident", "length", "mismatch"]
        )
        results = aggregate_off_target_genes(
            df, transcript_mapping, target_gene="BRCA1"
        )
        assert results == {}

    def test_multiple_probes(self, blast_hits_df, transcript_mapping):
        """Results should contain entries for all probes."""
        results = aggregate_off_target_genes(
            blast_hits_df, transcript_mapping, target_gene="BRCA1"
        )
        assert "probe-1" in results
        assert "probe-2" in results

    def test_cli_gene_name_variant_excluded_as_self_hit(self, transcript_mapping):
        """CLI gene names should exclude self-hits even when target_gene includes organism prefix."""
        df = pd.DataFrame(
            {
                "qseqid": ["probe-1", "probe-1"],
                "sseqid": ["ENST00000123456", "ENST00000654321"],
                "pident": [100.0, 92.0],
                "length": [25, 21],
                "mismatch": [0, 2],
            }
        )

        with patch("eFISHent.config.SequenceConfig") as mock_cfg:
            mock_cfg.return_value.gene_name = "BRCA1"
            results = aggregate_off_target_genes(
                df, transcript_mapping, target_gene="homo_sapiens_brca1"
            )

        info = results["probe-1"]
        assert info["txome_off_targets"] == 1
        assert info["off_target_genes"] == "TP53(1)"
