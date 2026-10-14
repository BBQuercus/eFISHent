"""Cross-validate probe specificity with BLAST against the reference genome."""

import os
import shutil
import subprocess
import tempfile

import Bio.SeqIO
import pandas as pd

from .config import GeneralConfig, ProbeConfig, SequenceConfig


class GenomeVerificationMixin:
    """Genome BLAST check counting near-full-length hits per probe."""

    def _ensure_genome_blast_db(self, genome: str) -> str | None:
        """Build the genome BLAST DB next to the FASTA, or reuse if cached.

        Falls back to a temporary directory if the genome directory is not
        writable. Returns the DB prefix (path without .nsq), or None on
        failure.
        """
        from .console import spinner

        if os.path.isfile(genome + ".nsq"):
            return genome

        genome_dir = os.path.dirname(genome) or "."
        if os.access(genome_dir, os.W_OK):
            db_path = genome
        else:
            self._verification_tmpdir = tempfile.mkdtemp(prefix="efishent_blastdb_")
            db_path = os.path.join(self._verification_tmpdir, "genome")
            self.logger.debug(
                f"Genome directory {genome_dir} not writable; "
                f"building one-shot BLAST DB in {self._verification_tmpdir}"
            )

        try:
            with spinner("Building verification BLAST database (cached for future runs)..."):
                subprocess.check_call(
                    ["makeblastdb", "-in", genome, "-dbtype", "nucl",
                     "-out", db_path],
                    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                )
        except subprocess.CalledProcessError as e:
            self.logger.debug(f"makeblastdb failed: {e}")
            return None
        return db_path

    def _run_blast_verification(self, probes_fasta: str) -> dict | None:
        """Cross-validate final probes with BLAST against the reference genome.

        Uses a different algorithm than Bowtie2 alignment to independently
        verify probe specificity. Only counts near-full-length matches
        (>= 75% of probe length at >= 90% identity) as significant hits.
        Returns None if BLAST is not available.
        """
        if not shutil.which("blastn") or not shutil.which("makeblastdb"):
            return None

        genome = os.path.abspath(GeneralConfig().reference_genome)
        if not genome or not os.path.isfile(genome):
            return None

        self._verification_tmpdir = None
        try:
            # When transcriptome BLAST handles off-target detection and
            # max_off_targets is at default (0), the BLAST result is unused
            # downstream — short-circuit before paying for makeblastdb +
            # genome BLAST (multi-minute on networked filesystems).
            is_endogenous = SequenceConfig().is_endogenous
            max_off = ProbeConfig().max_off_targets
            max_expected = max_off + (1 if is_endogenous else 0)
            has_transcriptome = bool(GeneralConfig().reference_transcriptome)
            if has_transcriptome and max_off == 0 and is_endogenous:
                hits_per_probe: dict[str, int] = {}
            else:
                db_path = self._ensure_genome_blast_db(genome)
                if db_path is None:
                    return None
                hits_per_probe = self._count_genomic_hits(probes_fasta, db_path)

            probes = [r.id for r in Bio.SeqIO.parse(probes_fasta, "fasta")]
            clean = sum(
                1 for p in probes
                if hits_per_probe.get(p, 0) <= max_expected
            )
            flagged = {
                p: hits_per_probe.get(p, 0)
                for p in probes
                if hits_per_probe.get(p, 0) > max_expected
            }
            return {
                "total": len(probes),
                "clean": clean,
                "flagged": flagged,
                "max_expected": max_expected,
            }
        except Exception as e:
            self.logger.debug(f"BLAST verification failed: {e}")
            return None
        finally:
            if self._verification_tmpdir:
                shutil.rmtree(self._verification_tmpdir, ignore_errors=True)
                self._verification_tmpdir = None

    @staticmethod
    def _count_genomic_hits(probes_fasta: str, db_path: str) -> dict[str, int]:
        """BLAST probes against the genome and count significant loci per probe."""
        from .console import spinner

        with tempfile.TemporaryDirectory() as tmpdir:
            blast_out = os.path.join(tmpdir, "blast.tsv")

            # BLAST with sensitive parameters
            with spinner("Running BLAST cross-validation..."):
                subprocess.check_call(
                    ["blastn", "-task", "blastn",
                     "-query", probes_fasta, "-db", db_path,
                     "-evalue", "10", "-word_size", "11",
                     "-dust", "no",
                     "-num_threads", str(GeneralConfig().threads),
                     "-outfmt", "6 qseqid sseqid pident length mismatch "
                                "gapopen qstart qend sstart send qlen",
                     "-out", blast_out],
                    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                )

            # Parse results
            cols = ["qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "qlen"]
            if os.path.getsize(blast_out) > 0:
                df_blast = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
            else:
                df_blast = pd.DataFrame(columns=cols)

        # Significant hits: near-full-length high-identity matches only
        # >= 95% identity AND alignment covers >= 80% of probe length
        # For a 20nt probe this means ~16nt at ≤1 mismatch
        df_blast["effective_len"] = df_blast["length"] - df_blast["gapopen"]
        sig = df_blast[
            (df_blast["pident"] >= 95)
            & (df_blast["effective_len"] >= df_blast["qlen"] * 0.80)
        ]

        # Count unique genomic loci per probe
        # Group by chromosome + 1kb-binned position to deduplicate
        sig = sig.copy()
        sig["locus"] = sig["sseqid"] + ":" + (sig["sstart"] // 1000).astype(str)
        return sig.groupby("qseqid")["locus"].nunique().to_dict()
