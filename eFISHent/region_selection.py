"""Select exonic or intronic regions of a gene sequence using the annotation."""

import logging
import os

import Bio.Seq
import Bio.SeqRecord

from .config import GeneralConfig


class RegionSelectionMixin:
    """Extract exon/intron regions of the target gene; mixed into PrepareSequence."""

    logger = logging.getLogger("custom-logger")

    def _extract_exon_boundaries(
        self, gene_name: str, annotation_path: str
    ) -> list[tuple[int, int]] | None:
        """Extract exon boundaries for a gene from GTF annotation.

        Returns a list of (start, end) tuples in transcript-relative coordinates,
        sorted by start position. Returns None if no exons found.
        """
        try:
            from .indexing import PrepareAnnotationFile, read_annotation
            parq_path = annotation_path + ".parq"
            if not os.path.isfile(parq_path):
                task = PrepareAnnotationFile()
                task.prepare_gtf_file(annotation_path, parq_path)
            df = read_annotation(
                parq_path,
                ["feature", "gene_name", "gene_id", "transcript_id", "start", "end"],
            )
        except Exception as e:
            self.logger.debug(f"Could not parse GTF for region selection: {e}")
            return None

        # Filter for exons of this gene
        gene_exons = df[
            (df["feature"] == "exon")
            & (df["gene_name"].str.upper() == gene_name.upper())
        ]

        if gene_exons.empty:
            # Try gene_id match
            gene_exons = df[
                (df["feature"] == "exon")
                & (df["gene_id"].str.contains(gene_name, case=False, na=False))
            ]

        if gene_exons.empty:
            self.logger.debug(f"No exons found for gene {gene_name} in GTF")
            return None

        # Get the transcript with the most exons (likely longest/primary)
        if "transcript_id" in gene_exons.columns:
            transcript_counts = gene_exons["transcript_id"].value_counts()
            primary_transcript = transcript_counts.index[0]
            gene_exons = gene_exons[gene_exons["transcript_id"] == primary_transcript]

        # Get gene start for relative coordinates
        gene_start = int(gene_exons["start"].min())

        exons = []
        for _, row in gene_exons.iterrows():
            start = int(row["start"]) - gene_start
            end = int(row["end"]) - gene_start
            exons.append((start, end))

        exons.sort()
        self.logger.debug(
            f"Found {len(exons)} exons for {gene_name}, "
            f"gene span: {gene_start}-{gene_start + exons[-1][1]}"
        )
        return exons

    def _apply_region_selection(
        self,
        sequence: Bio.SeqRecord.SeqRecord,
        target_regions: str,
        gene_name: str,
    ) -> Bio.SeqRecord.SeqRecord:
        """Extract specific regions (exon/intron) from a sequence.

        For 'exon' mode (default): returns the sequence as-is (transcriptome
        sequences are already exon-only).
        For 'intron' mode: requires genome + GTF to extract intronic regions.
        For 'both' mode: requires genome + GTF for the full pre-mRNA.
        """
        if target_regions == "exon":
            return sequence  # Default: transcript sequence is exon-only

        exons = self._region_selection_exons(target_regions, gene_name)
        if exons is None:
            return sequence

        if target_regions == "intron":
            return self._extract_intronic_sequence(sequence, exons, gene_name)

        if target_regions == "both":
            # Full pre-mRNA: the sequence is already the full gene
            self.logger.info(
                f"Using full pre-mRNA sequence ({len(str(sequence.seq)):,} nt)"
            )
        return sequence

    def _region_selection_exons(
        self, target_regions: str, gene_name: str
    ) -> list[tuple[int, int]] | None:
        """Get exon boundaries needed for region selection.

        Returns None (with a warning) if the required genome/annotation is
        missing or no exons were found, in which case selection falls back
        to exon-only.
        """
        annotation = GeneralConfig().reference_annotation
        genome = GeneralConfig().reference_genome

        if not annotation or not genome:
            self.logger.warning(
                f"--target-regions {target_regions} requires --reference-annotation "
                f"and --reference-genome. Falling back to exon-only."
            )
            return None

        exons = self._extract_exon_boundaries(gene_name, annotation)
        if exons is None:
            self.logger.warning(
                f"Could not find exon boundaries for {gene_name}. "
                f"Falling back to exon-only."
            )
        return exons

    def _extract_intronic_sequence(
        self,
        sequence: Bio.SeqRecord.SeqRecord,
        exons: list[tuple[int, int]],
        gene_name: str,
    ) -> Bio.SeqRecord.SeqRecord:
        """Concatenate the intronic regions (gaps between exons) of a sequence.

        Returns the sequence unchanged if there are no gaps between exons.
        """
        seq_str = str(sequence.seq)
        intronic_parts = []
        prev_end = 0
        for exon_start, exon_end in exons:
            if exon_start > prev_end:
                intronic_parts.append(seq_str[prev_end:exon_start])
            prev_end = exon_end

        if not intronic_parts:
            self.logger.warning(
                f"No intronic regions found for {gene_name}. "
                f"Gene may be single-exon."
            )
            return sequence

        intronic_seq = "".join(intronic_parts)
        self.logger.info(
            f"Extracted {len(intronic_parts)} intronic regions "
            f"({len(intronic_seq):,} nt) from {len(seq_str):,} nt gene"
        )
        return Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(intronic_seq),
            id=f"{sequence.id}_introns",
            name=sequence.name,
            description=f"intronic regions of {sequence.id}",
        )
