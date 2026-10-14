"""Report transcriptome off-targets and their expression risk for final probes."""

import os
import re
import shutil
import subprocess
import tempfile

import pandas as pd

from . import util
from .config import GeneralConfig, ProbeConfig
from .gene_annotation import aggregate_off_target_genes, build_transcript_gene_map

# Matches each "GENE(count)" entry of the off_target_genes column
_OFF_TARGET_GENE_RE = re.compile(r"([^,\s]+)\(\d+\)")


def _parse_off_target_genes(genes_str: str) -> list[str]:
    """Extract gene names (in order, with repeats) from an off_target_genes entry."""
    return [match.group(1) for match in _OFF_TARGET_GENE_RE.finditer(genes_str)]


def _iter_probe_off_target_genes(df: pd.DataFrame):
    """Yield (index, gene names) for each probe with off-target genes."""
    for idx, genes_str in df["off_target_genes"].items():
        if genes_str:
            yield idx, _parse_off_target_genes(genes_str)


def _classify_expression_risk(gene: str, gene_expression: dict[str, float]) -> str | None:
    """Classify expression risk for a single off-target gene."""
    if gene not in gene_expression:
        return None
    expr = gene_expression[gene]
    if expr > 50:
        return f"{gene}:HIGH({expr:.0f})"
    if expr > 10:
        return f"{gene}:moderate({expr:.0f})"
    if expr > 1:
        return f"{gene}:low({expr:.0f})"
    return None  # <1 TPM = negligible


class OffTargetReportMixin:
    """Transcriptome off-target details and expression risk columns."""

    def _add_transcriptome_details(self, df: pd.DataFrame) -> None:
        """Add transcriptome off-target details to probe table.

        Re-BLASTs final probes against the transcriptome and maps hits
        to gene names via GTF annotation.
        """
        txome = GeneralConfig().reference_transcriptome
        annotation = GeneralConfig().reference_annotation

        # Initialize empty columns
        df["txome_off_targets"] = 0
        df["off_target_genes"] = ""
        df["worst_match"] = ""

        if not txome:
            return

        # Check for existing txome hits CSV from filtering stage
        hits_csv = os.path.join(
            util.get_output_dir(),
            f"{util.get_gene_name()}_txome_hits.csv",
        )

        if not os.path.isfile(hits_csv):
            # Re-BLAST final probes against transcriptome
            hits_csv = self._blast_final_probes(df, txome)
            if not hits_csv:
                return

        try:
            df_hits = pd.read_csv(hits_csv)
        except Exception as e:
            self.logger.debug(f"Could not read txome hits: {e}")
            return

        if df_hits.empty:
            return

        # Apply filtering thresholds matching transcriptome_filter.py
        cfg = ProbeConfig()
        min_match_len = cfg.min_blast_match_length
        if min_match_len == 0:
            min_match_len = max(18, int(0.8 * cfg.min_length))

        if "effective_len" not in df_hits.columns:
            df_hits["effective_len"] = df_hits["length"] - df_hits["gapopen"]

        df_filtered = df_hits[
            (df_hits["effective_len"] >= min_match_len)
            & (df_hits["pident"] >= cfg.blast_identity_threshold)
        ]

        if df_filtered.empty:
            return

        # Build gene name mapping
        mapping: dict[str, str] = {}
        if annotation:
            # build_transcript_gene_map handles parquet discovery internally
            mapping = build_transcript_gene_map(annotation)

        target_gene = util.get_gene_name(hashed=False)
        results = aggregate_off_target_genes(df_filtered, mapping, target_gene)

        # Map results back to probe dataframe
        for idx, row in df.iterrows():
            probe_id = row["name"]
            if probe_id in results:
                info = results[probe_id]
                df.at[idx, "txome_off_targets"] = info["txome_off_targets"]
                df.at[idx, "off_target_genes"] = info["off_target_genes"]
                df.at[idx, "worst_match"] = info["worst_match"]

    def _blast_final_probes(
        self, df: pd.DataFrame, txome: str
    ) -> str | None:
        """BLAST final probes against the transcriptome for detailed reporting."""
        if not shutil.which("blastn"):
            return None

        txome_db = os.path.abspath(txome)
        # Check if BLAST DB exists
        if not os.path.isfile(txome_db + ".nsq"):
            return None

        from .console import spinner
        from .transcriptome_filter import (
            BLAST_EVALUE,
            BLAST_MAX_TARGETS,
            antisense_hits,
        )

        try:
            # Write final probes to temp FASTA
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            ) as tmp:
                for _, row in df.iterrows():
                    tmp.write(f">{row['name']}\n{row['sequence']}\n")
                tmp_fasta = tmp.name

            blast_out = tmp_fasta.replace(".fasta", "_blast.tsv")
            cfg = ProbeConfig()
            perc_identity = max(
                cfg.blast_identity_threshold,
                100.0 * 15 / cfg.max_length,
            )

            blast_task = "blastn-short" if cfg.max_length <= 30 else "blastn"
            args = [
                "blastn", "-task", blast_task,
                "-query", tmp_fasta,
                "-db", txome_db,
                "-evalue", BLAST_EVALUE,
                "-word_size", "7",
                "-gapopen", "5", "-gapextend", "2",
                "-reward", "1", "-penalty", "-3",
                "-dust", "no",
                "-max_target_seqs", BLAST_MAX_TARGETS,
                "-perc_identity", str(perc_identity),
                "-num_threads", str(GeneralConfig().threads),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                           "qstart qend sstart send evalue bitscore",
                "-out", blast_out,
            ]

            with spinner("BLASTing final probes against transcriptome..."):
                subprocess.check_call(
                    args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
                )

            # Parse and save as CSV
            columns = [
                "qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore",
            ]
            hits_csv = os.path.join(
                util.get_output_dir(),
                f"{util.get_gene_name()}_txome_hits.csv",
            )

            if os.path.getsize(blast_out) > 0:
                df_blast = pd.read_csv(
                    blast_out, sep="\t", header=None, names=columns
                )
                antisense_hits(df_blast).to_csv(hits_csv, index=False)
            else:
                pd.DataFrame(columns=columns).to_csv(hits_csv, index=False)

            # Clean up temp files
            os.remove(tmp_fasta)
            os.remove(blast_out)

            return hits_csv

        except Exception as e:
            self.logger.debug(f"Final probe BLAST failed: {e}")
            return None

    def _annotate_expression_risk(self, df: pd.DataFrame) -> None:
        """Add expression risk annotation based on off-target gene expression.

        When an expression table and GTF are provided, looks up expression
        levels for off-target genes and assigns risk categories.
        """
        df["expression_risk"] = ""

        annotation = GeneralConfig().reference_annotation
        count_table = ProbeConfig().encode_count_table

        if not annotation or "off_target_genes" not in df.columns:
            return

        gene_expression = self._load_gene_expression(count_table, annotation)

        # Annotate each probe
        for idx, genes in _iter_probe_off_target_genes(df):
            risks = [
                risk
                for risk in (
                    _classify_expression_risk(gene, gene_expression) for gene in genes
                )
                if risk
            ]
            if risks:
                df.at[idx, "expression_risk"] = "; ".join(risks)

    def _load_gene_expression(
        self, count_table: str, annotation: str
    ) -> dict[str, float]:
        """Build a gene->expression mapping from the count table, if available.

        Expression values are keyed by both gene ID and gene name (via GTF).
        """
        gene_expression: dict[str, float] = {}
        if not count_table or not os.path.isfile(count_table):
            return gene_expression

        try:
            sep = "\t" if count_table.endswith(".tsv") else ","
            df_expr = pd.read_csv(count_table, sep=sep)
            if len(df_expr.columns) < 2:
                return gene_expression

            # First col = gene_id, second col = expression value
            # Build gene_id -> gene_name mapping from GTF
            mapping = build_transcript_gene_map(annotation)

            # Map gene_ids to expression values, using gene_name as key
            for _, row in df_expr.iterrows():
                gid = str(row.iloc[0])
                expr_val = float(row.iloc[1])
                # Try to get gene name from mapping
                gname = mapping.get(gid, gid)
                gene_expression[gname] = expr_val
                gene_expression[gid] = expr_val
        except Exception as e:
            self.logger.debug(f"Could not load expression table: {e}")
        return gene_expression
