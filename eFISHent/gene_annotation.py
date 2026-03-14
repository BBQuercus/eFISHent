"""Map transcript IDs to gene names using GTF annotation files."""

import logging
import os
import re
from typing import Dict, Optional

import pandas as pd

logger = logging.getLogger("custom-logger")

# Cache for transcript-to-gene mapping
_transcript_gene_cache: Dict[str, Dict[str, str]] = {}


def build_transcript_gene_map(gtf_path: str) -> Dict[str, str]:
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

    mapping: Dict[str, str] = {}

    # Try parquet first (prepared by indexing.py)
    parquet_path = gtf_path
    if not gtf_path.endswith(".parquet"):
        parquet_path = os.path.splitext(gtf_path)[0] + ".parquet"

    if os.path.isfile(parquet_path):
        mapping = _parse_parquet_gtf(parquet_path)
    elif os.path.isfile(gtf_path):
        mapping = _parse_raw_gtf(gtf_path)
    else:
        logger.warning(f"GTF file not found: {gtf_path}")

    _transcript_gene_cache[gtf_path] = mapping
    logger.debug(f"Built transcript->gene mapping with {len(mapping)} entries")
    return mapping


def _parse_parquet_gtf(parquet_path: str) -> Dict[str, str]:
    """Parse a parquet-formatted GTF file for transcript->gene mapping."""
    df = pd.read_parquet(parquet_path)
    mapping: Dict[str, str] = {}

    # Check which columns are available
    has_gene_name = "gene_name" in df.columns
    has_gene_id = "gene_id" in df.columns
    has_transcript_id = "transcript_id" in df.columns

    if has_transcript_id and (has_gene_name or has_gene_id):
        name_col = "gene_name" if has_gene_name else "gene_id"
        pairs = df[["transcript_id", name_col]].dropna().drop_duplicates()
        for _, row in pairs.iterrows():
            tid = str(row["transcript_id"])
            gname = str(row[name_col])
            if tid and gname:
                mapping[tid] = gname
                # Also map versioned IDs (ENST00000123456.1 -> gene)
                base_id = tid.split(".")[0]
                if base_id != tid:
                    mapping[base_id] = gname

    # Also build gene_id -> gene_name mapping
    if has_gene_name and has_gene_id:
        gene_pairs = df[["gene_id", "gene_name"]].dropna().drop_duplicates()
        for _, row in gene_pairs.iterrows():
            gid = str(row["gene_id"])
            gname = str(row["gene_name"])
            if gid and gname and gid != gname:
                mapping[gid] = gname

    return mapping


def _parse_raw_gtf(gtf_path: str) -> Dict[str, str]:
    """Parse a raw GTF file for transcript->gene mapping."""
    import gzip

    mapping: Dict[str, str] = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:  # type: ignore
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            attrs = fields[8]

            # Extract transcript_id
            tid_match = re.search(r'transcript_id "([^"]+)"', attrs)
            if not tid_match:
                continue
            tid = tid_match.group(1)

            # Extract gene_name (preferred) or gene_id
            gname_match = re.search(r'gene_name "([^"]+)"', attrs)
            if gname_match:
                gname = gname_match.group(1)
            else:
                gid_match = re.search(r'gene_id "([^"]+)"', attrs)
                if gid_match:
                    gname = gid_match.group(1)
                else:
                    continue

            mapping[tid] = gname
            base_id = tid.split(".")[0]
            if base_id != tid:
                mapping[base_id] = gname

    return mapping


def map_transcript_to_gene(
    transcript_id: str, mapping: Dict[str, str]
) -> str:
    """Map a single transcript ID to its gene name.

    Handles various ID formats:
    - Ensembl: ENST00000123456.1
    - RefSeq: NM_001234.5
    - gffread: gene_id|transcript_id or similar pipe-delimited formats

    Args:
        transcript_id: The transcript/subject ID from BLAST output.
        mapping: Pre-built transcript->gene mapping dict.

    Returns:
        Gene name if found, otherwise the original transcript_id.
    """
    # Direct lookup
    if transcript_id in mapping:
        return mapping[transcript_id]

    # Try without version suffix
    base_id = transcript_id.split(".")[0]
    if base_id in mapping:
        return mapping[base_id]

    # Handle pipe-delimited IDs (gffread format: gene|transcript|...)
    if "|" in transcript_id:
        for part in transcript_id.split("|"):
            part = part.strip()
            if part in mapping:
                return mapping[part]
            part_base = part.split(".")[0]
            if part_base in mapping:
                return mapping[part_base]

    # Handle underscore-delimited IDs
    if "_" in transcript_id:
        parts = transcript_id.split("_")
        # Try first two parts (common in RefSeq: NM_001234)
        if len(parts) >= 2:
            refseq_id = f"{parts[0]}_{parts[1]}"
            if refseq_id in mapping:
                return mapping[refseq_id]

    return transcript_id


def aggregate_off_target_genes(
    blast_hits: pd.DataFrame,
    mapping: Dict[str, str],
    target_gene: str,
) -> Dict[str, Dict]:
    """Aggregate BLAST hits per probe into off-target gene summaries.

    Args:
        blast_hits: DataFrame with BLAST output columns (qseqid, sseqid, pident, length, etc.)
        mapping: transcript->gene mapping dict
        target_gene: Name of the target gene to exclude from off-targets

    Returns:
        Dict mapping probe_id -> {
            "txome_off_targets": int (unique off-target transcript count),
            "off_target_genes": str (comma-separated "GENE(count)" list),
            "worst_match": str (best match quality "pident%/length_bp/mismatches_mm"),
        }
    """
    target_lower = target_gene.lower()
    results: Dict[str, Dict] = {}

    for probe_id, group in blast_hits.groupby("qseqid"):
        gene_counts: Dict[str, int] = {}
        worst_pident = 0.0
        worst_length = 0
        worst_mm = 999

        for _, hit in group.iterrows():
            gene = map_transcript_to_gene(str(hit["sseqid"]), mapping)

            # Skip self-hits
            if gene.lower() == target_lower:
                continue
            # Also skip if the sseqid contains the target gene name
            if target_lower in str(hit["sseqid"]).lower():
                continue

            gene_counts[gene] = gene_counts.get(gene, 0) + 1

            # Track worst (highest quality) off-target match
            pident = float(hit["pident"])
            length = int(hit["length"])
            mm = int(hit.get("mismatch", 0))
            if pident > worst_pident or (pident == worst_pident and length > worst_length):
                worst_pident = pident
                worst_length = length
                worst_mm = mm

        # Format gene list
        if gene_counts:
            sorted_genes = sorted(gene_counts.items(), key=lambda x: -x[1])
            gene_str = ", ".join(f"{g}({c})" for g, c in sorted_genes)
            worst_str = f"{worst_pident:.0f}%/{worst_length}bp/{worst_mm}mm"
        else:
            gene_str = ""
            worst_str = ""

        results[str(probe_id)] = {
            "txome_off_targets": sum(gene_counts.values()),
            "off_target_genes": gene_str,
            "worst_match": worst_str,
        }

    return results
