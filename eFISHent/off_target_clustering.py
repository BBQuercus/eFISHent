"""Limit and warn about final probes clustering on the same off-target gene."""

import pandas as pd

from . import util
from .config import ProbeConfig
from .off_target_report import _iter_probe_off_target_genes


class OffTargetClusteringMixin:
    """Per-gene off-target cap and clustering/undruggable-target warnings."""

    def _apply_off_target_cap(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove probes to cap how many hit the same off-target gene.

        If max_probes_per_off_target is set and > 0, iteratively removes
        the lowest-quality probe hitting over-represented off-target genes.
        """
        cap = ProbeConfig().max_probes_per_off_target
        if cap <= 0 or "off_target_genes" not in df.columns:
            return df

        df = df.copy()
        max_iterations = len(df)  # safety limit

        for _ in range(max_iterations):
            # Find genes exceeding the cap
            over_cap = {
                gene: probes
                for gene, probes in self._probes_per_off_target_gene(df).items()
                if len(probes) > cap
            }

            if not over_cap:
                break

            # Remove the lowest-quality probe among those hitting the most over-represented gene
            worst_gene = max(over_cap, key=lambda g, oc=over_cap: len(oc[g]))
            probe_indices = over_cap[worst_gene]
            # Find the probe with lowest quality among these
            worst_idx = df.loc[probe_indices, "quality"].idxmin()
            self.logger.debug(
                f"Off-target cap: removing probe at idx {worst_idx} "
                f"(gene {worst_gene} hit by {len(probe_indices)} probes, cap={cap})"
            )
            df = df.drop(worst_idx)

        return df.reset_index(drop=True)

    @staticmethod
    def _probes_per_off_target_gene(df: pd.DataFrame) -> dict[str, list[int]]:
        """Map each off-target gene to the indices of the probes hitting it."""
        gene_probes: dict[str, list[int]] = {}
        for idx, genes in _iter_probe_off_target_genes(df):
            for gene in genes:
                gene_probes.setdefault(gene, []).append(idx)
        return gene_probes

    # Common gene family prefixes for paralog detection
    _GENE_FAMILY_PREFIXES = [
        "ACT", "MYH", "MYL", "HIST", "HLA", "KRT", "COL", "TUBB", "TUBA",
        "HSP", "RPL", "RPS", "EIF", "PCDH", "CDH", "CLDN",
    ]

    def _detect_off_target_clustering(self, df: pd.DataFrame) -> list[str]:
        """Detect when multiple probes hit the same off-target gene.

        Returns a list of warning messages. Also detects paralog families
        where the target may be inherently undruggable for FISH.
        """
        if "off_target_genes" not in df.columns:
            return []

        # Count how many probes hit each off-target gene
        gene_probe_counts: dict[str, int] = {}
        for _, genes in _iter_probe_off_target_genes(df):
            for gene in dict.fromkeys(genes):
                gene_probe_counts[gene] = gene_probe_counts.get(gene, 0) + 1

        cluster_threshold = 5
        n_probes = len(df)
        target_gene = util.get_gene_name(hashed=False)

        return [
            self._clustering_warning(gene, count, n_probes, target_gene)
            for gene, count in sorted(gene_probe_counts.items(), key=lambda x: -x[1])
            if count >= cluster_threshold
        ]

    @classmethod
    def _is_paralog(cls, gene: str, target_gene: str) -> bool:
        """Check whether two genes share a known gene family prefix."""
        return any(
            gene.upper().startswith(prefix) and target_gene.upper().startswith(prefix)
            for prefix in cls._GENE_FAMILY_PREFIXES
        )

    @classmethod
    def _clustering_warning(
        cls, gene: str, count: int, n_probes: int, target_gene: str
    ) -> str:
        """Build the warning for ``count`` probes clustering on off-target ``gene``."""
        is_paralog = cls._is_paralog(gene, target_gene)
        # Undruggable detection: ≥50% of probes hit same off-target
        is_undruggable = n_probes > 0 and count >= n_probes * 0.5

        if is_undruggable and is_paralog:
            return (
                f"UNDRUGGABLE TARGET: {count}/{n_probes} probes "
                f"({count * 100 // n_probes}%) hit {gene} "
                f"(paralog of {target_gene}). "
                f"This gene cannot be uniquely targeted by exonic FISH probes "
                f"due to conserved gene family homology. "
                f"Alternatives: (1) use --target-regions intron for intronic probes "
                f"(introns diverge between paralogs), "
                f"(2) choose a different target gene."
            )
        if is_undruggable:
            return (
                f"UNDRUGGABLE TARGET: {count}/{n_probes} probes "
                f"({count * 100 // n_probes}%) hit off-target gene {gene}. "
                f"This target may not be uniquely targetable. "
                f"Consider intronic probes (--target-regions intron) or "
                f"a different target gene."
            )
        if is_paralog:
            return (
                f"Off-target clustering: {count}/{n_probes} probes hit {gene} "
                f"(paralog of {target_gene}). This target may not be uniquely "
                f"targetable by FISH due to gene family homology. "
                f"Consider intronic probes instead."
            )
        return (
            f"Off-target clustering: {count}/{n_probes} probes hit "
            f"off-target gene {gene}."
        )
