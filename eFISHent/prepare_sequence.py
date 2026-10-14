"""Download Entrez Gene Probes from NCBI.

Select the right strand and select intronic/exonic regions.
"""

import logging
import os
import subprocess

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .config import GeneralConfig, SequenceConfig


class DownloadEntrezGeneSequence(luigi.Task):
    """Download genomic sequence as fasta from NCBI."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        fname = f"{util.get_gene_name()}_entrez.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def get_nuccore_query(
        self, ensembl_id: str, gene_name: str, organism_name: str
    ) -> str:
        """Build an NCBI nuccore query that returns only RefSeq transcripts.

        Searches the nuccore database directly with RefSeq and biomol filters
        to avoid returning genomic sequences, ESTs, or other non-transcript records.
        For Ensembl IDs, first resolves to an NCBI Gene ID via the gene database.
        """
        refseq_filter = "refseq[filter] AND (biomol_mrna[prop] OR biomol_rna[prop])"

        if ensembl_id:
            # Resolve Ensembl ID to NCBI Gene ID via the gene database
            ensembl_query = f"({ensembl_id})"
            if organism_name:
                ensembl_query += f" AND {organism_name}[Organism]"

            search = subprocess.run(
                ["esearch", "-db", "gene", "-query", ensembl_query],
                check=True,
                capture_output=True,
            )
            uid_result = subprocess.run(
                ["efetch", "-format", "uid"],
                input=search.stdout,
                check=True,
                capture_output=True,
            )
            gene_id = uid_result.stdout.decode().strip()
            if not gene_id:
                raise LookupError(
                    f"No NCBI Gene found for Ensembl ID '{ensembl_id}'. "
                    "Please check the ID and try again."
                )
            # Take first Gene ID if multiple returned
            gene_id = gene_id.split("\n")[0].strip()
            self.logger.debug(f"Resolved Ensembl ID '{ensembl_id}' to Gene ID {gene_id}.")
            return f"{gene_id}[Gene ID] AND {refseq_filter}"

        if gene_name and organism_name:
            return (
                f"({gene_name}[Gene Name]) AND {organism_name}[Organism]"
                f" AND {refseq_filter}"
            )

        raise ValueError(
            "For downloading Entrez Gene Probes, "
            "you need to specify the gene name and organism name or provide an Ensembl ID."
        )

    def fetch_entrez(
        self, ensembl_id: str, gene_name: str, organism_name: str
    ) -> str:
        """Retrieve the fasta sequence data from NCBI.

        Queries the nuccore database directly with RefSeq transcript filters
        to ensure only curated mRNA/ncRNA sequences are returned.
        """
        query = self.get_nuccore_query(ensembl_id, gene_name, organism_name)
        self.logger.debug(f"Fetching from Entrez nuccore using query '{query}'.")

        search = subprocess.run(
            ["esearch", "-db", "nuccore", "-query", query],
            check=True,
            capture_output=True,
        )
        fetch = subprocess.run(
            ["efetch", "-format", "fasta"],
            input=search.stdout,
            check=True,
            capture_output=True,
        )
        fasta = fetch.stdout.decode()

        # Check for empty results or missing FASTA content
        if not fasta or ">" not in fasta:
            raise LookupError(
                "No RefSeq transcript found for your query. "
                "Please check your gene name and organism and try again. "
                "Or provide a custom fasta file with --sequence-file."
            )
        return fasta

    def run(self):
        util.log_stage_start(self.logger, "DownloadEntrezGeneSequence")
        if not SequenceConfig().is_plus_strand:
            self.logger.warning(

                    "Downloading from entrez but using negative strand as target! "
                    "Entrez usually provides sequences 5'-3'. "
                    "Proceed with caution."

            )
        fasta = self.fetch_entrez(
            SequenceConfig().ensembl_id,
            SequenceConfig().gene_name,
            SequenceConfig().organism_name,
        )

        # Log which sequences were retrieved
        headers = [line for line in fasta.split("\n") if line.startswith(">")]
        n_records = len(headers)
        if n_records == 1:
            self.logger.info(f"Downloaded 1 transcript: {headers[0][1:]}")
        else:
            self.logger.info(f"Downloaded {n_records} transcript isoforms from NCBI")
            max_show = 3
            for header in headers[:max_show]:
                self.logger.debug(f"  {header[1:]}")
            if n_records > max_show:
                self.logger.debug(f"  ...and {n_records - max_show} more")

        with self.output().open("w") as outfile:
            outfile.write(fasta)


class PrepareSequence(luigi.Task):
    """Prepare gene sequence.

    Select right strand and exon/intron/both.
    """

    logger = logging.getLogger("custom-logger")

    def requires(self):
        tasks = {}
        if not SequenceConfig().sequence_file:
            tasks["entrez"] = DownloadEntrezGeneSequence()
        return tasks

    def output(self):
        fname = f"{util.get_gene_name()}_sequence.fasta"
        return luigi.LocalTarget(os.path.join(util.get_output_dir(), fname))

    def select_sequence(
        self, sequences: list[Bio.SeqRecord.SeqRecord]
    ) -> Bio.SeqRecord.SeqRecord:
        """Select a single sequence as template.

        When multiple records are present (e.g. transcript isoforms),
        filters to only keep those matching the target gene name
        (NCBI sometimes returns related genes), then selects the longest.
        """
        import re

        if not sequences:
            raise ValueError(
                "No records found in fasta file. "
                "Please ensure at least one sequence is present."
            )

        # Filter by target gene name to exclude unrelated genes
        # (e.g. NCBI returns POTEF when searching for ACTB)
        gene_name = SequenceConfig().gene_name
        if gene_name and len(sequences) > 1:
            pattern = re.compile(
                rf'\b{re.escape(gene_name)}\b', re.IGNORECASE
            )
            matching = [
                s for s in sequences
                if pattern.search(s.description)
            ]
            if matching and len(matching) < len(sequences):
                self.logger.info(
                    f"Filtered {len(sequences) - len(matching)} sequences "
                    f"not matching gene name '{gene_name}'."
                )
                sequences = matching

        if len(sequences) > 1:
            longest = max(sequences, key=lambda s: len(s.seq))
            self.logger.info(
                f"Found {len(sequences)} transcript isoforms. "
                f"Selected longest: {longest.id} ({len(longest.seq):,} nt)."
            )
            return longest
        return sequences[0]

    def select_strand(
        self, sequence: Bio.SeqRecord.SeqRecord, is_plus_strand: bool
    ) -> Bio.SeqRecord.SeqRecord:
        """Select right strand for probes to target."""
        if is_plus_strand:
            sequence = sequence.reverse_complement(id=True, description=True)
            self.logger.debug("Converted sequence to reverse complement.")
        return sequence

    @staticmethod
    def is_fasta_formatted(fname):
        if not os.path.exists(fname):
            return False

        with open(fname) as fasta_file:
            in_sequence = False
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    in_sequence = True
                elif not in_sequence or not all(
                    base in "ATCGNatcgn" for base in line
                ):
                    return False
            return True

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

    def run(self):
        util.log_stage_start(self.logger, "PrepareSequence")
        input_file = (
            self.input()["entrez"].path
            if "entrez" in self.input()
            else SequenceConfig().sequence_file
        )
        if not self.is_fasta_formatted(input_file):
            raise ValueError(
                "Fasta file is incorrectly formatted. "
                "Check wikipedia: https://en.wikipedia.org/wiki/FASTA_format"
            )
        sequences = list(Bio.SeqIO.parse(input_file, format="fasta"))
        sequence = self.select_sequence(sequences)
        sequence = self.select_strand(sequence, SequenceConfig().is_plus_strand)

        # Apply region selection if specified
        target_regions = getattr(SequenceConfig(), "target_regions", "exon")
        if target_regions != "exon":
            gene_name = SequenceConfig().gene_name or util.get_gene_name(hashed=False)
            sequence = self._apply_region_selection(
                sequence, target_regions, gene_name
            )

        self.logger.info(f'Selected sequence and strand called "{sequence.id}".')
        Bio.SeqIO.write(sequence, self.output().path, format="fasta")
