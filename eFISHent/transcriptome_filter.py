"""Filter probes based on BLAST hits against a reference transcriptome.

BLAST parameters adapted from TrueProbes (Neuert Lab, 2025). Uses blastn
(not megablast) with word_size 7 to detect short partial matches on probes.
The approach of BLASTing probes against the transcriptome (rather than genome)
catches off-target binding to expressed transcripts that alignment-only methods miss.
"""

import logging
import math
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import Bio.SeqIO
import luigi
import pandas as pd

from . import util
from .alignment import AlignProbeCandidates
from .config import GeneralConfig
from .config import ProbeConfig


def _run_blast_chunk(args):
    """Run a single BLAST process (top-level function for pickling)."""
    subprocess.check_call(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


class BuildTranscriptomeBlastDB(luigi.Task):
    """Build a BLAST nucleotide database from a transcriptome FASTA."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        txome = os.path.abspath(GeneralConfig().reference_transcriptome)
        return luigi.LocalTarget(txome + ".nsq")

    def run(self):
        from .console import spinner

        txome = os.path.abspath(GeneralConfig().reference_transcriptome)
        args = [
            "makeblastdb",
            "-in", txome,
            "-dbtype", "nucl",
            "-out", txome,
        ]
        self.logger.debug(f"Running makeblastdb with - {' '.join(args)}")
        with spinner("Building BLAST database..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )
        self.logger.info("Finished building BLAST transcriptome database.")


class TranscriptomeFiltering(luigi.Task):
    """Filter probes with significant off-target hits in the transcriptome."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return {
            "probes": AlignProbeCandidates(),
            "blastdb": BuildTranscriptomeBlastDB(),
        }

    def output(self):
        return {
            "fasta": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(),
                    f"{util.get_gene_name()}_txome.fasta",
                )
            ),
            "table": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(),
                    f"{util.get_gene_name()}_txome_hits.csv",
                )
            ),
        }

    @staticmethod
    def _blast_args(task, query, db, perc_identity, threads, out):
        """Build BLAST command-line arguments."""
        return [
            "blastn",
            "-task", task,
            "-query", query,
            "-db", db,
            "-evalue", "100",
            "-word_size", "7",
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "1",
            "-penalty", "-3",
            "-dust", "no",
            "-max_target_seqs", "50",
            "-perc_identity", str(perc_identity),
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore",
            "-out", out,
        ]

    def run(self):
        util.log_stage_start(self.logger, "TranscriptomeFiltering")
        probe_fasta = self.input()["probes"]["fasta"].path
        sequences = list(Bio.SeqIO.parse(probe_fasta, "fasta"))

        txome_db = os.path.abspath(GeneralConfig().reference_transcriptome)
        blast_out = os.path.splitext(probe_fasta)[0] + "_blast.tsv"

        config = ProbeConfig()
        max_probe_len = config.max_length
        perc_identity = max(config.blast_identity_threshold, 100.0 * 15 / max_probe_len)

        # Use blastn-short for probe-length queries (<30bp) — optimized seed
        # strategy and scoring for short sequences
        blast_task = "blastn-short" if max_probe_len <= 30 else "blastn"
        threads = GeneralConfig().threads

        from .console import spinner

        # Parallel BLAST: split probes into chunks and run independent BLAST
        # processes. This scales better than BLAST's own -num_threads because
        # each process gets its own I/O and memory pipeline.
        n_chunks = max(1, min(threads, math.ceil(len(sequences) / 500)))

        if n_chunks <= 1:
            # Small probe set — single BLAST call
            args = self._blast_args(
                blast_task, probe_fasta, txome_db, perc_identity, threads, blast_out,
            )
            self.logger.debug(f"Running BLAST with - {' '.join(args)}")
            with spinner("BLASTing probes against transcriptome..."):
                subprocess.check_call(
                    args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
                )
        else:
            # Split probes into chunks and BLAST in parallel
            self.logger.debug(
                f"Splitting {len(sequences)} probes into {n_chunks} chunks "
                f"for parallel BLAST"
            )
            chunk_size = math.ceil(len(sequences) / n_chunks)
            tmpdir = tempfile.mkdtemp(prefix="efishent_blast_")
            chunk_files = []
            chunk_outs = []

            try:
                for i in range(n_chunks):
                    chunk = sequences[i * chunk_size : (i + 1) * chunk_size]
                    if not chunk:
                        continue
                    chunk_fa = os.path.join(tmpdir, f"chunk_{i}.fasta")
                    chunk_out = os.path.join(tmpdir, f"chunk_{i}.tsv")
                    Bio.SeqIO.write(chunk, chunk_fa, "fasta")
                    chunk_files.append(chunk_fa)
                    chunk_outs.append(chunk_out)

                # Each chunk gets 1 thread — parallelism comes from processes
                all_args = [
                    self._blast_args(
                        blast_task, cf, txome_db, perc_identity, 1, co,
                    )
                    for cf, co in zip(chunk_files, chunk_outs)
                ]

                with spinner(
                    f"BLASTing probes against transcriptome ({n_chunks} parallel)..."
                ):
                    with ProcessPoolExecutor(max_workers=n_chunks) as pool:
                        list(pool.map(_run_blast_chunk, all_args))

                # Concatenate results
                with open(blast_out, "w") as fout:
                    for co in chunk_outs:
                        if os.path.isfile(co) and os.path.getsize(co) > 0:
                            with open(co) as fin:
                                fout.write(fin.read())
            finally:
                import shutil
                shutil.rmtree(tmpdir, ignore_errors=True)

        # Parse BLAST output
        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore",
        ]

        if os.path.getsize(blast_out) > 0:
            df_blast = pd.read_csv(blast_out, sep="\t", header=None, names=columns)
        else:
            df_blast = pd.DataFrame(columns=columns)

        # Post-BLAST filtering (per TrueProbes approach):
        # effective match length = alignment length - gaps
        # Default min_match_len: max(18, 0.8 * min_probe_length)
        min_match_len = config.min_blast_match_length
        if min_match_len == 0:
            min_match_len = max(18, int(0.8 * config.min_length))
        self.logger.debug(f"Using min_blast_match_length = {min_match_len}")
        df_blast["effective_len"] = df_blast["length"] - df_blast["gapopen"]
        df_hits = df_blast[
            (df_blast["effective_len"] >= min_match_len)
            & (df_blast["pident"] >= config.blast_identity_threshold)
        ]

        # Count off-target hits per probe (exclude self-hits)
        # Strategy 1: gene name in subject ID (works for gffread-style headers)
        # Strategy 2: GTF-based transcript-to-gene mapping (works for Ensembl IDs)
        gene_name_raw = util.get_gene_name(hashed=False).lower()
        gene_name = re.escape(gene_name_raw)
        gene_pattern = rf"(?:^|[|_\-.\s]){gene_name}(?:$|[|_\-.\s])"

        # Also try extracting a cleaner gene name by stripping common suffixes
        # (e.g., "EIF2B1_cds" -> "EIF2B1", "METTL3_cds_odd" -> "METTL3")
        gene_name_clean = re.sub(r"[_-]?(cds|mrna|transcript|seq|odd|even).*$", "", gene_name_raw)
        if gene_name_clean != gene_name_raw:
            gene_name_clean_esc = re.escape(gene_name_clean)
            gene_pattern_clean = rf"(?:^|[|_\-.\s]){gene_name_clean_esc}(?:$|[|_\-.\s])"
        else:
            gene_pattern_clean = None

        # Include the raw --gene-name CLI value (e.g., "URA3") which may differ
        # from the file-based gene_name_raw (e.g., "saccharomyces_cerevisiae_ura3")
        from .config import SequenceConfig
        cli_gene_name = (SequenceConfig().gene_name or "").strip().lower()

        # Build transcript-to-gene mapping from GTF if available
        self_transcript_ids = set()
        gtf_path = GeneralConfig().reference_annotation
        if gtf_path:
            try:
                from .gene_annotation import build_transcript_gene_map
                tx_map = build_transcript_gene_map(gtf_path)
                # Find all transcript IDs that belong to the target gene
                # Try raw name, cleaned name, and the CLI --gene-name value
                search_names = {gene_name_raw, gene_name_clean}
                if cli_gene_name:
                    search_names.add(cli_gene_name)
                for search_name in search_names:
                    for tid, gname in tx_map.items():
                        if gname.lower() == search_name:
                            self_transcript_ids.add(tid.lower())
                if self_transcript_ids:
                    self.logger.debug(
                        f"GTF mapping: found {len(self_transcript_ids)} transcripts "
                        f"for target gene '{gene_name_clean}'"
                    )
            except Exception as e:
                self.logger.debug(f"GTF mapping failed, using name-based exclusion: {e}")

        # Build regex for CLI gene name (e.g., "URA3" -> matches "YEL021W_URA3_mRNA")
        if cli_gene_name:
            cli_gene_esc = re.escape(cli_gene_name)
            cli_gene_pattern = rf"(?:^|[|_\-.\s]){cli_gene_esc}(?:$|[|_\-.\s])"
        else:
            cli_gene_pattern = None

        def is_self_hit(sseqid):
            """Check if a subject sequence ID is a self-hit."""
            sid_lower = sseqid.lower()
            base = sid_lower.split(".")[0]

            # Check gene name patterns (gffread-style, cleaned, CLI)
            patterns = [p for p in (gene_pattern, gene_pattern_clean, cli_gene_pattern) if p]
            return (
                any(re.search(p, sid_lower) for p in patterns)
                or sid_lower in self_transcript_ids
                or base in self_transcript_ids
            )

        off_target_counts = {}
        for probe_id, group in df_hits.groupby("qseqid"):
            off_targets = group[~group["sseqid"].apply(is_self_hit)]
            off_target_counts[probe_id] = len(off_targets["sseqid"].unique())

        # Filter probes
        max_off = config.max_transcriptome_off_targets
        candidates = [
            seq for seq in sequences
            if off_target_counts.get(seq.id, 0) <= max_off
        ]

        # Cross-hybridization check: probes with long contiguous match at
        # >=95% identity to off-target transcripts. These near-perfect matches
        # to abundant or spatially concentrated RNAs create noisy images.
        # Threshold: match must cover >=85% of the probe length (per-probe).
        # For a 20nt probe: >=17nt. For a 24nt probe: >=20nt.
        cross_hyb_min_frac = 0.85
        cross_hyb_override = config.cross_hyb_min_length
        # Build probe length lookup from sequences
        probe_lengths = {seq.id: len(seq.seq) for seq in sequences}
        df_cross_hyb_base = df_blast[
            (df_blast["pident"] >= 95.0)
            & (df_blast["gapopen"] == 0)  # contiguous = no gaps
        ]
        cross_hyb_probes = set()
        if not df_cross_hyb_base.empty:
            for probe_id, group in df_cross_hyb_base.groupby("qseqid"):
                probe_len = probe_lengths.get(probe_id, config.min_length)
                if cross_hyb_override > 0:
                    threshold = cross_hyb_override
                else:
                    threshold = max(15, int(cross_hyb_min_frac * probe_len))
                hits = group[group["effective_len"] >= threshold]
                if hits.empty:
                    continue
                off_targets = hits[~hits["sseqid"].apply(is_self_hit)]
                if not off_targets.empty and off_targets["sseqid"].nunique() > 0:
                    cross_hyb_probes.add(probe_id)

        if cross_hyb_probes:
            if config.reject_cross_hybridization:
                self.logger.debug(
                    f"Rejecting {len(cross_hyb_probes)} probes with strong "
                    f"cross-hybridization (>=85% of probe length contiguous at >=95% identity)"
                )
                candidates = [
                    seq for seq in candidates if seq.id not in cross_hyb_probes
                ]
            else:
                self.logger.debug(
                    f"Cross-hybridization warning: {len(cross_hyb_probes)} probes have "
                    f">=85% contiguous match at >=95% identity to off-target transcripts"
                )

        # Save hits table
        df_hits.to_csv(self.output()["table"].path, index=False)

        util.log_and_check_candidates(
            self.logger, "TranscriptomeFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, format="fasta")

        # Clean up temp blast output
        if os.path.exists(blast_out):
            os.remove(blast_out)
