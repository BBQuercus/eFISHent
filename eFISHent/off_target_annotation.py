"""Off-target filters for alignments based on the genome annotation."""

import bisect
import logging
import os

import pandas as pd

from . import constants
from .config import ProbeConfig
from .indexing import read_annotation


class AnnotationOffTargetFilterMixin:
    """Filter alignment hits by annotated gene loci and expression.

    Mixed into AlignProbeCandidates, which provides the task inputs and
    filter settings.
    """

    logger = logging.getLogger("custom-logger")
    encode_count_table: str

    @property
    def count_table(self) -> pd.DataFrame:
        """Read and verify a RNAseq normalized count table."""
        extension = os.path.splitext(self.fname_count)[1].lower()
        if extension not in (".tsv", ".csv", ".txt"):
            raise ValueError(
                "The count table must be provided as TSV, CSV, or TXT file. "
                f"Only found - {self.fname_count}"
            )
        sep = "\t" if extension in (".tsv", ".txt") else ","
        df = pd.read_csv(self.fname_count, sep=sep, index_col=0).reset_index().dropna()
        if len(df.columns) < 2:
            raise ValueError(
                "The count table must contain at least 2 columns."
                f"Only found - {df.columns}"
            )
        self.logger.debug(f"Mapping original the first two columns of - {df.columns}.")
        self.logger.debug(f"Read count table with {len(df)} entries.")
        self.logger.debug(f"Found dtypes - {df.dtypes}")
        self.logger.debug(f"First column - {df.iloc[0]}")

        df = df.iloc[:, :2]
        df.columns = ["gene_id", "count"]
        try:
            df = df.astype({"gene_id": str, "count": float})
        except ValueError as e:
            raise ValueError(
                "Could not convert count table columns to the right type. "
                "Please ensure the first two colums are of type string and number/float/int, respectively. "
            ) from e

        # Remove potential version numbers
        df["gene_id"] = df["gene_id"].apply(lambda x: x.split(".")[0])
        return df

    @property
    def gtf_table(self) -> pd.DataFrame:
        """Parse prepared parquet-formatted gtf file."""
        return read_annotation(self.fname_gtf, constants.GTF_COLUMNS)[constants.GTF_COLUMNS]

    @property
    def norm_table(self) -> pd.DataFrame:
        """Merge gtf and count table."""
        df = pd.merge(
            self.count_table, self.gtf_table, how="left", on="gene_id",
            validate="many_to_many",
        ).dropna()
        return df

    def get_most_expressed_genes(self, percentage: float) -> pd.Series:
        """Get gene IDs expressed above percentile threshold."""
        threshold = self.norm_table["count"].quantile(1 - percentage / 100)
        return self.norm_table[self.norm_table["count"] >= threshold]["gene_id"]

    def _gene_loci_with_biotype(
        self, biotypes: set[str], filter_name: str, loci_name: str
    ) -> pd.DataFrame | None:
        """Get unique gene loci (seqname, start, end) annotated with any of biotypes.

        Returns None (and logs why) if the annotation has no biotype column or
        no gene matches, in which case the calling filter should be skipped.
        """
        self.fname_gtf = self.input()["annotation"].path
        gtf_full = read_annotation(
            self.fname_gtf,
            ["seqname", "start", "end", "feature", "gene_biotype", "gene_type"],
        )

        # Find the biotype column (Ensembl: gene_biotype, GENCODE: gene_type)
        biotype_col = next(
            (col for col in ["gene_biotype", "gene_type"] if col in gtf_full.columns),
            None,
        )
        if biotype_col is None:
            self.logger.debug(
                f"No gene_biotype/gene_type column in GTF — skipping {filter_name} filter"
            )
            return None

        loci = gtf_full[
            (gtf_full[biotype_col].isin(biotypes))
            & (gtf_full["feature"] == "gene")
        ][["seqname", "start", "end"]].drop_duplicates()

        if loci.empty:
            self.logger.debug(f"No {loci_name} found in GTF annotation")
            return None
        return loci

    @staticmethod
    def _sorted_intervals_by_chrom(loci: pd.DataFrame) -> dict:
        """Build a bisect-ready interval lookup from gene loci.

        Returns dict of {chrom: (sorted_starts, ends)} with intervals sorted by
        (start, end).
        """
        by_chrom = {}
        for _, gene in loci.iterrows():
            chrom = str(gene["seqname"])
            if chrom not in by_chrom:
                by_chrom[chrom] = ([], [])
            by_chrom[chrom][0].append(int(gene["start"]))
            by_chrom[chrom][1].append(int(gene["end"]))

        for chrom, (starts, ends) in by_chrom.items():
            paired = sorted(zip(starts, ends))
            by_chrom[chrom] = (
                [s for s, _ in paired],
                [e for _, e in paired],
            )
        return by_chrom

    def filter_rrna_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove probes with off-target hits on ribosomal RNA genes.

        rRNA is ~80% of cellular RNA — even weak binding causes intense
        background. Checks gene_biotype/gene_type for rRNA and Mt_rRNA.
        Any probe with at least one rRNA hit is vetoed entirely.
        """
        rrna_types = {"rRNA", "Mt_rRNA", "rRNA_pseudogene"}
        rrna_genes = self._gene_loci_with_biotype(
            rrna_types, filter_name="rRNA", loci_name="rRNA genes"
        )
        if rrna_genes is None:
            return df_sam

        self.logger.debug(f"Found {len(rrna_genes)} rRNA gene regions in GTF")

        # Build sorted interval lookup per chromosome for efficient overlap check
        rrna_by_chrom = {}
        for _, gene in rrna_genes.iterrows():
            chrom = gene["seqname"]
            if chrom not in rrna_by_chrom:
                rrna_by_chrom[chrom] = ([], [])
            rrna_by_chrom[chrom][0].append(gene["start"])
            rrna_by_chrom[chrom][1].append(gene["end"])

        def _hits_rrna(row):
            chrom = row["rname"]
            if chrom not in rrna_by_chrom:
                return False
            pos = int(row["pos"])
            starts, ends = rrna_by_chrom[chrom]
            return any(s <= pos <= e for s, e in zip(starts, ends))

        hits_rrna = df_sam.apply(_hits_rrna, axis=1)
        probes_to_remove = set(df_sam.loc[hits_rrna, "qname"])

        if probes_to_remove:
            self.logger.debug(
                f"Removing {len(probes_to_remove)} probes hitting rRNA off-targets"
            )
            df_sam = df_sam[~df_sam["qname"].isin(probes_to_remove)]

        return df_sam

    def filter_intergenic_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove off-target hits in intergenic regions before counting.

        Hits not overlapping any annotated gene are unlikely to produce RNA
        that the FISH probe could bind, so they are excluded from off-target
        counts. Uses the prepared GTF annotation.
        """
        self.fname_gtf = self.input()["annotation"].path
        gtf = self.gtf_table
        genes = gtf[gtf["feature"] == "gene"]

        total_before = len(df_sam)

        def _overlaps_gene(row):
            chrom = row["rname"]
            pos = int(row["pos"])
            return (
                (genes["seqname"] == chrom)
                & (genes["start"] <= pos)
                & (genes["end"] >= pos)
            ).any()

        is_genic = df_sam.apply(_overlaps_gene, axis=1)
        df_sam = df_sam[is_genic]
        removed = total_before - len(df_sam)
        self.logger.debug(
            f"Removed {removed} alignments in intergenic regions"
        )
        return df_sam

    def filter_pseudogene_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove off-target hits overlapping annotated pseudogene loci.

        Pseudogenes (especially processed pseudogenes from retrotransposition)
        produce little to no RNA, so alignment hits to these regions are not
        meaningful off-targets for RNA-FISH probes. Removing them prevents
        pseudogene-rich genes (e.g. ACTB, GAPDH) from losing all probes.
        """
        pseudogene_types = {
            "processed_pseudogene", "unprocessed_pseudogene", "pseudogene",
            "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
            "translated_processed_pseudogene", "polymorphic_pseudogene",
            "IG_pseudogene", "TR_pseudogene",
            "rRNA_pseudogene", "unitary_pseudogene",
            "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
            "TR_J_pseudogene", "TR_V_pseudogene",
        }
        pseudogenes = self._gene_loci_with_biotype(
            pseudogene_types, filter_name="pseudogene", loci_name="pseudogene regions"
        )
        if pseudogenes is None:
            return df_sam

        self.logger.debug(
            f"Found {len(pseudogenes)} pseudogene regions in GTF"
        )

        pseudo_by_chrom = self._sorted_intervals_by_chrom(pseudogenes)

        def _hits_pseudogene(row):
            chrom = row["rname"]
            if chrom not in pseudo_by_chrom:
                return False
            pos = int(row["pos"])
            starts, ends = pseudo_by_chrom[chrom]
            idx = bisect.bisect_right(starts, pos) - 1
            return any(
                starts[i] <= pos <= ends[i]
                for i in range(max(0, idx), min(len(starts), idx + 2))
            )

        total_before = len(df_sam)
        is_pseudo = df_sam.apply(_hits_pseudogene, axis=1)
        df_sam = df_sam[~is_pseudo]
        removed = total_before - len(df_sam)
        self.logger.debug(
            f"Removed {removed} alignments overlapping pseudogene loci"
        )
        return df_sam

    def filter_expression_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove probes aligning to highly expressed off-target genes."""
        self.fname_count = self.encode_count_table
        self.fname_gtf = self.input()["annotation"].path

        expressed_genes = self.get_most_expressed_genes(
            ProbeConfig().max_expression_percentage
        )
        self.logger.debug(
            f"Found {len(expressed_genes)} highly expressed genes "
            f"(top {ProbeConfig().max_expression_percentage}%)"
        )

        # Join alignments with gene annotations by position overlap
        gtf = self.gtf_table
        gtf = gtf[gtf["feature"] == "gene"]

        # For each alignment, find overlapping genes
        probes_to_remove = set()
        for _, row in df_sam.iterrows():
            chrom = row["rname"]
            pos = int(row["pos"])
            # Find genes overlapping this position
            overlapping = gtf[
                (gtf["seqname"] == chrom)
                & (gtf["start"] <= pos)
                & (gtf["end"] >= pos)
            ]
            for _, gene in overlapping.iterrows():
                gene_id = gene["gene_id"]
                if gene_id in expressed_genes.values:
                    probes_to_remove.add(row["qname"])

        if probes_to_remove:
            self.logger.debug(
                f"Removing {len(probes_to_remove)} probes hitting highly expressed off-targets"
            )
            df_sam = df_sam[~df_sam["qname"].isin(probes_to_remove)]

        return df_sam
