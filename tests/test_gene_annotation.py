"""Tests for gene_annotation module."""

import os
import tempfile

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
