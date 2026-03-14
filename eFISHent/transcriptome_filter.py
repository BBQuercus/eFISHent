"""Filter probes based on BLAST hits against a reference transcriptome.

BLAST parameters adapted from TrueProbes (Neuert Lab, 2025). Uses blastn
(not megablast) with word_size 7 to detect short partial matches on probes.
The approach of BLASTing probes against the transcriptome (rather than genome)
catches off-target binding to expressed transcripts that alignment-only methods miss.
"""

import logging
import os
import re
import subprocess

import Bio.SeqIO
import luigi
import pandas as pd

from . import util
from .alignment import AlignProbeCandidates
from .config import GeneralConfig
from .config import ProbeConfig


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

    def run(self):
        util.log_stage_start(self.logger, "TranscriptomeFiltering")
        probe_fasta = self.input()["probes"]["fasta"].path
        sequences = list(Bio.SeqIO.parse(probe_fasta, "fasta"))

        txome_db = os.path.abspath(GeneralConfig().reference_transcriptome)
        blast_out = os.path.splitext(probe_fasta)[0] + "_blast.tsv"

        config = ProbeConfig()
        max_probe_len = config.max_length
        perc_identity = max(config.blast_identity_threshold, 100.0 * 15 / max_probe_len)

        # BLAST parameters derived from TrueProbes
        args = [
            "blastn",
            "-task", "blastn",
            "-query", probe_fasta,
            "-db", txome_db,
            "-evalue", "1000",
            "-word_size", "7",
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "1",
            "-penalty", "-3",
            "-dust", "no",
            "-num_alignments", "1000",
            "-perc_identity", str(perc_identity),
            "-num_threads", str(GeneralConfig().threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-out", blast_out,
        ]
        self.logger.debug(f"Running BLAST with - {' '.join(args)}")
        from .console import spinner

        with spinner("BLASTing probes against transcriptome..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )

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

        # Count off-target hits per probe (exclude self-hits by checking
        # if the probe maps to its own target — we identify self-hits by
        # the gene name appearing as a word in the subject sequence ID)
        gene_name = re.escape(util.get_gene_name(hashed=False).lower())
        # Use word boundary or delimiter-boundary matching to avoid substring
        # false positives (e.g., gene "NAG" matching transcript "GNAG1")
        gene_pattern = rf"(?:^|[|_\-.\s]){gene_name}(?:$|[|_\-.\s])"
        off_target_counts = {}
        for probe_id, group in df_hits.groupby("qseqid"):
            # Filter out self-hits where sseqid contains the target gene name
            # as a delimited word (not as a substring)
            off_targets = group[
                ~group["sseqid"].str.lower().str.contains(
                    gene_pattern, na=False, regex=True
                )
            ]
            off_target_counts[probe_id] = len(off_targets["sseqid"].unique())

        # Filter probes
        max_off = config.max_transcriptome_off_targets
        candidates = [
            seq for seq in sequences
            if off_target_counts.get(seq.id, 0) <= max_off
        ]

        # Secondary check: flag probes with >=16nt contiguous match at >=95% identity
        # These are warnings, not hard rejections
        cross_hyb_threshold = 16
        df_cross_hyb = df_blast[
            (df_blast["effective_len"] >= cross_hyb_threshold)
            & (df_blast["pident"] >= 95.0)
            & (df_blast["gapopen"] == 0)  # contiguous = no gaps
        ]
        if not df_cross_hyb.empty:
            cross_hyb_probes = set()
            for probe_id, group in df_cross_hyb.groupby("qseqid"):
                off_targets = group[
                    ~group["sseqid"].str.lower().str.contains(
                        gene_pattern, na=False, regex=True
                    )
                ]
                if len(off_targets["sseqid"].unique()) > 0:
                    cross_hyb_probes.add(probe_id)
            if cross_hyb_probes:
                self.logger.debug(
                    f"Cross-hybridization warning: {len(cross_hyb_probes)} probes have "
                    f">=16nt contiguous match at >=95% identity to off-target transcripts"
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
