"""Map transcript IDs to gene names using GTF annotation files."""

import logging
import os
import re

import pandas as pd

from .off_target_genes import (
    aggregate_off_target_genes as aggregate_off_target_genes,
)
from .off_target_genes import (
    map_transcript_to_gene as map_transcript_to_gene,
)

logger = logging.getLogger("custom-logger")

# Cache for transcript-to-gene mapping
_transcript_gene_cache: dict[str, dict[str, str]] = {}


def build_transcript_gene_map(gtf_path: str) -> dict[str, str]:
    """Build a mapping from transcript_id to gene_name using a GTF file.

    Supports Ensembl GTF (gene_name attribute), GENCODE GTF, and
    gffread-generated FASTAs (gene_id in header). Caches results
    per GTF path to avoid re-parsing.

    Args:
        gtf_path: Path to GTF or prepared parquet file.

    Returns:
        Dictionary mapping transcript_id -> gene_name.
    """
    if gtf_path in _transcript_gene_cache:
        return _transcript_gene_cache[gtf_path]

    mapping: dict[str, str] = {}

    # Try parquet first (prepared by indexing.py / PrepareAnnotationFile)
    # PrepareAnnotationFile outputs "{gtf_path}.parq", so check both extensions
    parquet_candidates = []
    if gtf_path.endswith((".parquet", ".parq")):
        parquet_candidates.append(gtf_path)
    else:
        parquet_candidates.append(gtf_path + ".parq")  # PrepareAnnotationFile output
        parquet_candidates.append(os.path.splitext(gtf_path)[0] + ".parquet")

    parquet_path = None
    for candidate in parquet_candidates:
        if os.path.isfile(candidate):
            parquet_path = candidate
            break

    if parquet_path is not None:
        mapping = _parse_parquet_gtf(parquet_path)
    elif os.path.isfile(gtf_path):
        mapping = _parse_raw_gtf(gtf_path)
    else:
        logger.warning(f"GTF file not found: {gtf_path}")

    _transcript_gene_cache[gtf_path] = mapping
    logger.debug(f"Built transcript->gene mapping with {len(mapping)} entries")
    return mapping


def _add_transcript(mapping: dict[str, str], tid: str, gname: str) -> None:
    """Map a transcript ID and its unversioned base ID to a gene name.

    Versioned IDs (ENST00000123456.1 -> gene) are also mapped without version.
    """
    mapping[tid] = gname
    base_id = tid.split(".")[0]
    if base_id != tid:
        mapping[base_id] = gname


def _string_pairs(df: pd.DataFrame, key_col: str, value_col: str) -> list[tuple[str, str]]:
    """Unique non-null (key, value) pairs of two columns, as strings."""
    pairs = df[[key_col, value_col]].dropna().drop_duplicates()
    return [(str(row[key_col]), str(row[value_col])) for _, row in pairs.iterrows()]


def _parse_parquet_gtf(parquet_path: str) -> dict[str, str]:
    """Parse a parquet-formatted GTF file for transcript->gene mapping."""
    from .indexing import read_annotation

    df = read_annotation(parquet_path, ["transcript_id", "gene_id", "gene_name"])
    mapping: dict[str, str] = {}
    columns = set(df.columns)

    name_col = "gene_name" if "gene_name" in columns else "gene_id"
    if {"transcript_id", name_col} <= columns:
        for tid, gname in _string_pairs(df, "transcript_id", name_col):
            if tid and gname:
                _add_transcript(mapping, tid, gname)

    # Also build gene_id -> gene_name mapping
    if {"gene_id", "gene_name"} <= columns:
        for gid, gname in _string_pairs(df, "gene_id", "gene_name"):
            if gid and gname and gid != gname:
                mapping[gid] = gname

    return mapping


def _transcript_and_gene(attrs: str) -> tuple[str, str] | None:
    """Extract transcript_id and gene_name (preferred) or gene_id from GTF attributes."""
    tid_match = re.search(r'transcript_id "([^"]+)"', attrs)
    if not tid_match:
        return None
    gene_match = re.search(r'gene_name "([^"]+)"', attrs) or re.search(
        r'gene_id "([^"]+)"', attrs
    )
    if not gene_match:
        return None
    return tid_match.group(1), gene_match.group(1)


def _parse_raw_gtf(gtf_path: str) -> dict[str, str]:
    """Parse a raw GTF file for transcript->gene mapping."""
    import gzip

    mapping: dict[str, str] = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:  # type: ignore
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            pair = _transcript_and_gene(fields[8])
            if pair is not None:
                _add_transcript(mapping, *pair)

    return mapping
