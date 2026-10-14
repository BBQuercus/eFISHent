"""Summarize BLAST off-target hits per probe by the genes they hit."""

import logging
import re

import pandas as pd

logger = logging.getLogger("custom-logger")


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
