"""Map transcript IDs to gene names using GTF annotation files."""

import logging
import os
import re

import pandas as pd

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


def _candidate_ids(transcript_id: str):
    """Yield the IDs to look up for a BLAST subject ID, in order of preference."""
    # Direct lookup, then without version suffix
    yield transcript_id
    yield transcript_id.split(".")[0]

    # Handle pipe-delimited IDs (gffread format: gene|transcript|...)
    if "|" in transcript_id:
        for part in transcript_id.split("|"):
            part = part.strip()
            yield part
            yield part.split(".")[0]

    # Handle underscore-delimited IDs: try first two parts
    # (common in RefSeq: NM_001234)
    parts = transcript_id.split("_")
    if len(parts) >= 2:
        yield f"{parts[0]}_{parts[1]}"


def map_transcript_to_gene(
    transcript_id: str, mapping: dict[str, str]
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
    for candidate in _candidate_ids(transcript_id):
        if candidate in mapping:
            return mapping[candidate]
    return transcript_id


def _cli_gene_name() -> str:
    """The --gene-name CLI value (lowercase), or "" if unavailable."""
    try:
        from .config import SequenceConfig

        return (SequenceConfig().gene_name or "").strip().lower()
    except Exception as e:
        logger.debug(f"Could not read --gene-name for self-hit detection: {e}")
        return ""


def _target_name_variants(target_gene: str) -> set[str]:
    """Build the set of target gene name variants for self-hit detection.

    target_gene may be "homo_sapiens_actb" (file-based) while the GTF
    gene name is just "ACTB", so also include the CLI --gene-name value
    and any suffix of the target after underscores (organism prefix).
    """
    target_lower = target_gene.lower()
    # Strip common suffixes to get clean gene name for matching
    # e.g., "EIF2B1_cds" -> "eif2b1", "METTL3_cds_odd" -> "mettl3"
    target_clean = re.sub(
        r"[_\-]?(cds|mrna|transcript|seq|odd|even).*$", "", target_lower
    )
    target_names = {target_lower, target_clean}
    # Add the bare gene name (last component after organism prefix)
    # e.g., "homo_sapiens_actb" -> "actb"
    parts = target_lower.split("_")
    if len(parts) > 1:
        target_names.add(parts[-1])
    cli_name = _cli_gene_name()
    if cli_name:
        target_names.add(cli_name)
    target_names.discard("")
    return target_names


def _is_target_hit(gene: str, sseqid: str, target_names: set[str]) -> bool:
    """Check if a hit is to the target gene (by gene name or subject ID)."""
    if gene.lower() in target_names:
        return True
    sseqid_lower = sseqid.lower()
    return any(name in sseqid_lower for name in target_names)


def _summarize_probe_hits(
    group: pd.DataFrame, mapping: dict[str, str], target_names: set[str]
) -> dict:
    """Summarize one probe's off-target BLAST hits."""
    gene_counts: dict[str, int] = {}
    worst_pident = 0.0
    worst_length = 0
    worst_mm = 999

    for _, hit in group.iterrows():
        gene = map_transcript_to_gene(str(hit["sseqid"]), mapping)
        if _is_target_hit(gene, str(hit["sseqid"]), target_names):
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

    gene_str = ""
    worst_str = ""
    if gene_counts:
        sorted_genes = sorted(gene_counts.items(), key=lambda x: -x[1])
        gene_str = ", ".join(f"{g}({c})" for g, c in sorted_genes)
        worst_str = f"{worst_pident:.0f}%/{worst_length}bp/{worst_mm}mm"

    return {
        "txome_off_targets": sum(gene_counts.values()),
        "off_target_genes": gene_str,
        "worst_match": worst_str,
    }


def aggregate_off_target_genes(
    blast_hits: pd.DataFrame,
    mapping: dict[str, str],
    target_gene: str,
) -> dict[str, dict]:
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
    target_names = _target_name_variants(target_gene)
    return {
        str(probe_id): _summarize_probe_hits(group, mapping, target_names)
        for probe_id, group in blast_hits.groupby("qseqid")
    }
