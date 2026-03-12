"""Clean up the output directory.

Remove the files that are not needed and prettify the kept output.
"""

from typing import Dict, List, Optional
import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import luigi
import pandas as pd

from . import util
from .alignment import AlignProbeCandidates
from .basic_filtering import get_gc_content
from .basic_filtering import get_melting_temp
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .kmers import BuildJellyfishIndex
from .kmers import get_max_kmer_counts_batch
from .optimization import OptimizeProbeCoverage
from .secondary_structure import get_free_energy


class CleanUpOutput(luigi.Task):
    """Clean up the output files and remove the intermediaries that are not needed."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_blast_required = (
            ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
        )

    def requires(self):
        tasks = {
            "optimize": OptimizeProbeCoverage(),
            "jellyfish": BuildJellyfishIndex(),
            "alignment": AlignProbeCandidates(),
        }
        return tasks

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("fasta", f"{util.get_gene_name()}.fasta"),
                ("table", f"{util.get_gene_name()}.csv"),
                ("config", f"{util.get_gene_name()}.txt"),
            ]
        }

    def prettify_table(
        self,
        sequences: List[Bio.SeqRecord.SeqRecord],
        basename: str,
        jellyfish_path: str,
        alignment_path: str,
        config: luigi.Config = ProbeConfig,  # type: ignore
    ) -> pd.DataFrame:
        """Create table with probe information."""
        # Create basic table with probe start/end positions and name
        df = util.create_data_table(sequences)
        sequences = sorted(sequences, key=lambda x: int(x.id.split("-")[-1]))

        # Add data columns
        df["GC"] = [round(get_gc_content(seq.seq), 2) for seq in sequences]
        df["TM"] = [
            round(
                get_melting_temp(
                    seq.seq, config().na_concentration, config().formamide_concentration  # type: ignore
                ),
                2,
            )
            for seq in sequences
        ]
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            df["deltaG"] = pool.map(get_free_energy, sequences)
        df["kmers"] = get_max_kmer_counts_batch(sequences, jellyfish_path)
        df_counts = pd.read_csv(alignment_path)
        df_counts = df_counts.groupby("qname", as_index=False).size()
        df["count"] = pd.merge(
            df, df_counts, how="left", left_on="name", right_on="qname",
            validate="many_to_one",
        )["size"].fillna(0)
        df["count"] -= 1

        # Compute per-probe quality score (0-100)
        df["quality"] = self._compute_quality_scores(df)

        # Create new/clean names
        df["name"] = [f"{basename}-{idx + 1}" for idx in df.index]
        return df

    @staticmethod
    def _compute_quality_scores(df: pd.DataFrame) -> pd.Series:
        """Compute a composite quality score (0-100) for each probe.

        Scores are based on how close each metric is to its ideal value:
        - Tm: closeness to midpoint of min/max range (30%)
        - GC: closeness to 50% (20%)
        - deltaG: less negative is better, 0 is ideal (20%)
        - kmers: lower is better (15%)
        - off-target count: lower is better (15%)
        """
        cfg = ProbeConfig()
        scores = pd.DataFrame(index=df.index)

        # Tm score: distance from midpoint of allowed range
        tm_mid = (cfg.min_tm + cfg.max_tm) / 2
        tm_range = max((cfg.max_tm - cfg.min_tm) / 2, 1)
        scores["tm"] = 1.0 - (df["TM"] - tm_mid).abs().clip(upper=tm_range) / tm_range

        # GC score: distance from 50%
        scores["gc"] = 1.0 - (df["GC"] - 50.0).abs() / 50.0

        # deltaG score: 0 is ideal, more negative is worse
        # Normalize against config threshold
        dg_max = abs(cfg.max_deltag) if cfg.max_deltag < 0 else 10.0
        scores["dg"] = 1.0 - df["deltaG"].abs().clip(upper=dg_max) / dg_max

        # K-mer score: 0 is ideal, higher is worse
        kmer_max = max(cfg.max_kmers, 1)
        scores["kmer"] = 1.0 - df["kmers"].clip(upper=kmer_max) / kmer_max

        # Off-target score: 0 is ideal
        ot_max = max(cfg.max_off_targets + 1, 1)
        scores["ot"] = 1.0 - df["count"].clip(lower=0, upper=ot_max) / ot_max

        # Clip all component scores to [0, 1] for robustness
        for col in scores.columns:
            scores[col] = scores[col].clip(lower=0.0, upper=1.0)

        # Weighted composite
        quality = (
            scores["tm"] * 0.30
            + scores["gc"] * 0.20
            + scores["dg"] * 0.20
            + scores["kmer"] * 0.15
            + scores["ot"] * 0.15
        )
        return (quality * 100).round(1)

    def _get_gene_length(self) -> int:
        """Get the actual gene sequence length from the sequence file."""
        # Try the prepared sequence file (available before intermediate cleanup)
        sequence_fasta = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_sequence.fasta"
        )
        if os.path.isfile(sequence_fasta):
            seq = next(Bio.SeqIO.parse(sequence_fasta, "fasta"))
            return len(seq)
        # Try user-provided sequence file
        if SequenceConfig().sequence_file and os.path.isfile(SequenceConfig().sequence_file):
            seq = next(Bio.SeqIO.parse(SequenceConfig().sequence_file, "fasta"))
            return len(seq)
        return 0

    def _compute_summary(self, df: pd.DataFrame) -> dict:
        """Compute summary statistics from the final probe DataFrame."""
        from .console import get_funnel_data

        gene_name = util.get_gene_name(hashed=False)

        # Get actual gene length, fall back to probe span
        gene_length = self._get_gene_length()
        if gene_length == 0 and len(df) > 0:
            gene_length = int(df["end"].max())

        # Compute covered base pairs by merging overlapping intervals
        covered_bp = 0
        if len(df) > 0:
            intervals = sorted(zip(df["start"], df["end"]))
            merged_start, merged_end = intervals[0]
            for s, e in intervals[1:]:
                if s <= merged_end:
                    merged_end = max(merged_end, e)
                else:
                    covered_bp += merged_end - merged_start
                    merged_start, merged_end = s, e
            covered_bp += merged_end - merged_start

        coverage_pct = (covered_bp / gene_length * 100) if gene_length > 0 else 0.0

        # Get initial candidate count from funnel data
        funnel = get_funnel_data()
        initial_count = funnel[0][1] if funnel else None

        return {
            "gene_name": gene_name,
            "probe_count": len(df),
            "initial_count": initial_count,
            "coverage_pct": coverage_pct,
            "gene_length": gene_length,
            "length_range": (int(df["length"].min()), int(df["length"].max())),
            "length_median": int(df["length"].median()),
            "tm_range": (float(df["TM"].min()), float(df["TM"].max())),
            "tm_median": float(df["TM"].median()),
            "gc_range": (float(df["GC"].min()), float(df["GC"].max())),
            "gc_median": float(df["GC"].median()),
        }

    def _run_blast_verification(self, probes_fasta: str) -> Optional[Dict]:
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

        from .console import spinner

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                db_path = os.path.join(tmpdir, "genome")
                blast_out = os.path.join(tmpdir, "blast.tsv")

                # Build temporary BLAST DB
                with spinner("Building verification BLAST database..."):
                    subprocess.check_call(
                        ["makeblastdb", "-in", genome, "-dbtype", "nucl",
                         "-out", db_path],
                        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                    )

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
                hits_per_probe = sig.groupby("qseqid")["locus"].nunique().to_dict()

                probes = [r.id for r in Bio.SeqIO.parse(probes_fasta, "fasta")]

                # For exogenous probes, 0 self-hits expected (not in genome)
                # For endogenous probes, 1 self-hit expected (target locus)
                is_endogenous = SequenceConfig().is_endogenous
                max_expected = ProbeConfig().max_off_targets + (1 if is_endogenous else 0)

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

    def prettify_sequences(self, df: pd.DataFrame) -> List[Bio.SeqRecord.SeqRecord]:
        """Clean up sequence id's and descriptions."""
        return [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(seq), id=name, name=name, description=""
            )
            for seq, name in zip(df["sequence"], df["name"])
        ]

    def prettify_section(self, section: str) -> str:
        """Prettify a single configuration section."""
        pretty = [f"{section.replace('Config', '')} configuration:"]
        for key, value in self.config[section].items():
            if value and value != '""':
                pretty.append(f"  {key.replace('_', '-')}: {value}")
        return "\n".join(pretty)

    def prettify_configuration(self) -> str:
        """Create a pretty file with all set/default configuration."""
        self.config = luigi.configuration.get_config()
        pretty = [self.prettify_section(section) for section in self.config.sections()]
        command = " ".join(sys.argv[1:])
        pretty.append(f"Command:\n  efishent {command}")
        return "\n\n".join(pretty)

    def remove_intermediates(self) -> None:
        """Remove intermediary files no longer used."""
        intermediary_files = glob.glob(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_*")
        )
        self.logger.debug(f"Found file intermediates - {intermediary_files}")
        try:
            intermediary_files.remove(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_entrez.fasta"
                )
            )
        except ValueError:
            self.logger.debug("Didn't find entrez file...")
        for filename in intermediary_files:
            self.logger.debug(f"Removing {filename}")
            os.remove(filename)

    def run(self):
        util.log_stage_start(self.logger, "CleanUpOutput")
        sequences = list(
            Bio.SeqIO.parse(self.input()["optimize"]["probes"].path, "fasta")  # type: ignore
        )
        df = self.prettify_table(
            sequences,
            basename=util.get_gene_name(hashed=False),
            jellyfish_path=self.input()["jellyfish"].path,  # type: ignore
            alignment_path=self.input()["alignment"]["table"].path,  # type: ignore
        )
        sequences = self.prettify_sequences(df)
        config = self.prettify_configuration()

        util.log_and_check_candidates(self.logger, "CleanUpOutput", len(sequences))

        # Store summary stats and probe data for completion message
        self._summary = self._compute_summary(df)
        self._probe_df = df

        df.to_csv(self.output()["table"].path, index=False)
        Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")
        with open(self.output()["config"].path, "w") as f:
            f.write(config)
        self.logger.debug(f'Saving all files using hash "{util.get_gene_name()}"')

        # Run BLAST verification on final probes
        self._verification = self._run_blast_verification(
            self.output()["fasta"].path
        )

        # Files to be deleted
        if not RunConfig().save_intermediates:
            self.remove_intermediates()
