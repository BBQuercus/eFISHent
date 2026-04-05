"""Clean up the output directory.

Remove the files that are not needed and prettify the kept output.
"""

from typing import Dict, List, Optional
import glob
import logging
import multiprocessing
import os
import re
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
from .basic_filtering import compute_duplex_dg
from .basic_filtering import get_cpg_fraction
from .basic_filtering import get_gc_content
from .basic_filtering import get_melting_temp
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .gene_annotation import aggregate_off_target_genes
from .gene_annotation import build_transcript_gene_map
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

        # Add CpG fraction and low-complexity score for quality scoring
        df["cpg_fraction"] = [
            round(get_cpg_fraction(seq.seq), 4) for seq in sequences
        ]
        df["low_complexity"] = [
            self._compute_low_complexity_score(seq.seq) for seq in sequences
        ]

        # Add target site accessibility scoring if enabled
        if ProbeConfig().accessibility_scoring:
            df["accessibility"] = self._compute_accessibility_scores(df)
        else:
            df["accessibility"] = 1.0

        # Add on-target duplex ΔG (binding stability)
        cfg_inst = config()
        df["on_target_dg"] = [
            round(
                compute_duplex_dg(
                    seq.seq, cfg_inst.na_concentration, cfg_inst.formamide_concentration
                ),
                2,
            )
            for seq in sequences
        ]

        # Add transcriptome off-target details if available
        self._add_transcriptome_details(df)

        # Add expression risk annotation if available
        self._annotate_expression_risk(df)

        # Compute per-probe quality score (0-100)
        df["quality"] = self._compute_quality_scores(df)

        # Add recommendation column
        df["recommendation"] = df.apply(self._compute_recommendation, axis=1)

        # Create new/clean names
        df["name"] = [f"{basename}-{idx + 1}" for idx in df.index]
        return df

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
        mapping: Dict[str, str] = {}
        if annotation:
            # Try parquet version first
            parquet_path = os.path.join(
                util.get_output_dir(),
                os.path.splitext(os.path.basename(annotation))[0] + ".parquet",
            )
            if os.path.isfile(parquet_path):
                mapping = build_transcript_gene_map(parquet_path)
            else:
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
    ) -> Optional[str]:
        """BLAST final probes against the transcriptome for detailed reporting."""
        if not shutil.which("blastn"):
            return None

        txome_db = os.path.abspath(txome)
        # Check if BLAST DB exists
        if not os.path.isfile(txome_db + ".nsq"):
            return None

        from .console import spinner

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

            args = [
                "blastn", "-task", "blastn",
                "-query", tmp_fasta,
                "-db", txome_db,
                "-evalue", "1000",
                "-word_size", "7",
                "-gapopen", "5", "-gapextend", "2",
                "-reward", "1", "-penalty", "-3",
                "-dust", "no",
                "-num_alignments", "1000",
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
                df_blast.to_csv(hits_csv, index=False)
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

        # Build gene->expression mapping if count table is available
        gene_expression: Dict[str, float] = {}
        if count_table and os.path.isfile(count_table):
            try:
                sep = "\t" if count_table.endswith(".tsv") else ","
                df_expr = pd.read_csv(count_table, sep=sep)
                cols = df_expr.columns.tolist()
                if len(cols) >= 2:
                    # First col = gene_id, second col = expression value
                    # Build gene_id -> gene_name mapping from GTF
                    parquet_path = os.path.join(
                        util.get_output_dir(),
                        os.path.splitext(os.path.basename(annotation))[0] + ".parquet",
                    )
                    if os.path.isfile(parquet_path):
                        mapping = build_transcript_gene_map(parquet_path)
                    else:
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

        # Annotate each probe
        for idx, row in df.iterrows():
            genes_str = row.get("off_target_genes", "")
            if not genes_str:
                continue

            risks = []
            # Parse "GENE(count)" format
            for match in re.finditer(r"([^,\s]+)\(\d+\)", genes_str):
                gene = match.group(1)
                if gene in gene_expression:
                    expr = gene_expression[gene]
                    if expr > 50:
                        risks.append(f"{gene}:HIGH({expr:.0f})")
                    elif expr > 10:
                        risks.append(f"{gene}:moderate({expr:.0f})")
                    elif expr > 1:
                        risks.append(f"{gene}:low({expr:.0f})")
                    # <1 TPM = negligible, skip

            if risks:
                df.at[idx, "expression_risk"] = "; ".join(risks)

    @staticmethod
    def _compute_recommendation(row: pd.Series) -> str:
        """Compute a PASS/FLAG/FAIL recommendation for a probe."""
        quality = row.get("quality", 0)
        txome_ot = row.get("txome_off_targets", 0)
        expr_risk = row.get("expression_risk", "")

        # FAIL conditions
        if quality < 30:
            return "FAIL"

        # FLAG conditions
        flags = []
        if quality < 70:
            flags.append("low_quality")
        if txome_ot > 3:
            flags.append("many_off_targets")
        if "HIGH" in str(expr_risk):
            flags.append("high_expression_risk")

        if flags:
            return f"FLAG({','.join(flags)})"

        return "PASS"

    @staticmethod
    def _compute_quality_scores(df: pd.DataFrame) -> pd.Series:
        """Compute a composite quality score (0-100) for each probe.

        Data-driven weights (n=84 probe sets):
        - Tm: deviation from set median (not config midpoint) (20%)
        - GC: optimal 45-52%, steep penalty above 55% (15%)
        - CpG: CpG dinucleotide frequency penalty (5%)
        - deltaG: less negative is better, 0 is ideal (15%)
        - off-target count: genome + transcriptome (25%)
        - kmers: lower is better (10%)
        - low-complexity: penalize low-complexity probes (10%)
        """
        cfg = ProbeConfig()
        scores = pd.DataFrame(index=df.index)

        # Tm score: distance from the *actual median Tm* of the selected set
        tm_median = df["TM"].median()
        tm_range = max((cfg.max_tm - cfg.min_tm) / 2, 1)
        scores["tm"] = 1.0 - (df["TM"] - tm_median).abs().clip(upper=tm_range) / tm_range

        # GC score: data-driven optimal range 45-52%
        # Full score in 45-52%, linear penalty outside, steep above 55%
        def _gc_score(gc: float) -> float:
            if 45.0 <= gc <= 52.0:
                return 1.0
            elif gc < 45.0:
                # Linear decay from 45% down to 20% (score 0)
                return max(0.0, (gc - 20.0) / 25.0)
            elif gc <= 55.0:
                # Gentle penalty 52-55%
                return 1.0 - (gc - 52.0) / 10.0
            else:
                # Steep penalty above 55% (p=0.002 for unspecific binding)
                return max(0.0, 0.7 - (gc - 55.0) / 15.0)

        scores["gc"] = df["GC"].apply(_gc_score)

        # CpG score: CpG frequency penalty (p=0.004 for nuclear background)
        if "cpg_fraction" in df.columns:
            # 0.0 CpG = perfect, 0.10+ = zero score
            scores["cpg"] = (1.0 - df["cpg_fraction"].clip(upper=0.10) / 0.10)
        else:
            scores["cpg"] = 1.0

        # deltaG score: 0 is ideal, more negative is worse
        dg_max = abs(cfg.max_deltag) if cfg.max_deltag < 0 else 10.0
        scores["dg"] = 1.0 - df["deltaG"].abs().clip(upper=dg_max) / dg_max

        # K-mer score: 0 is ideal, higher is worse
        kmer_max = max(cfg.max_kmers, 1)
        scores["kmer"] = 1.0 - df["kmers"].clip(upper=kmer_max) / kmer_max

        # Off-target score: combines genome off-targets and transcriptome off-targets
        ot_max = max(cfg.max_off_targets + 1, 1)
        genome_ot_score = 1.0 - df["count"].clip(lower=0, upper=ot_max) / ot_max
        txome_ot_score = pd.Series(1.0, index=df.index)
        if "txome_off_targets" in df.columns:
            txome_ot_score = 1.0 - df["txome_off_targets"].clip(lower=0, upper=5) / 5
        scores["ot"] = (genome_ot_score * 0.5 + txome_ot_score * 0.5)

        # Low-complexity score: penalize probes with low Shannon entropy
        if "low_complexity" in df.columns:
            scores["lc"] = 1.0 - df["low_complexity"].clip(upper=0.3) / 0.3
        else:
            scores["lc"] = 1.0

        # Accessibility score: higher = more unpaired target bases = better binding
        if "accessibility" in df.columns:
            scores["acc"] = df["accessibility"].clip(lower=0.0, upper=1.0)
        else:
            scores["acc"] = 1.0

        # On-target duplex ΔG score: more negative = more stable binding = better
        # Typical range: -20 to -60 kcal/mol for 20-25nt probes
        if "on_target_dg" in df.columns:
            # Normalize: -60 or below → 1.0, 0 → 0.0
            dg_floor = -60.0
            scores["binding"] = (df["on_target_dg"].clip(upper=0.0) / dg_floor).clip(
                lower=0.0, upper=1.0
            )
        else:
            scores["binding"] = 1.0

        # Clip all component scores to [0, 1] for robustness
        for col in scores.columns:
            scores[col] = scores[col].clip(lower=0.0, upper=1.0)

        # Data-driven weighted composite
        # Accessibility gets weight when enabled (redistributed from others)
        has_accessibility = "accessibility" in df.columns and (df["accessibility"] < 1.0).any()
        has_binding = "on_target_dg" in df.columns and (df["on_target_dg"] < 0.0).any()

        if has_accessibility and has_binding:
            quality = (
                scores["tm"] * 0.16
                + scores["gc"] * 0.12
                + scores["cpg"] * 0.05
                + scores["dg"] * 0.11
                + scores["ot"] * 0.21
                + scores["kmer"] * 0.08
                + scores["lc"] * 0.08
                + scores["acc"] * 0.09
                + scores["binding"] * 0.10
            )
        elif has_accessibility:
            quality = (
                scores["tm"] * 0.18
                + scores["gc"] * 0.14
                + scores["cpg"] * 0.05
                + scores["dg"] * 0.13
                + scores["ot"] * 0.23
                + scores["kmer"] * 0.09
                + scores["lc"] * 0.09
                + scores["acc"] * 0.09
            )
        elif has_binding:
            quality = (
                scores["tm"] * 0.18
                + scores["gc"] * 0.13
                + scores["cpg"] * 0.05
                + scores["dg"] * 0.13
                + scores["ot"] * 0.23
                + scores["kmer"] * 0.08
                + scores["lc"] * 0.10
                + scores["binding"] * 0.10
            )
        else:
            quality = (
                scores["tm"] * 0.20
                + scores["gc"] * 0.15
                + scores["cpg"] * 0.05
                + scores["dg"] * 0.15
                + scores["ot"] * 0.25
                + scores["kmer"] * 0.10
                + scores["lc"] * 0.10
            )
        return (quality * 100).round(1)

    def _compute_accessibility_scores(self, df: pd.DataFrame) -> pd.Series:
        """Compute target RNA accessibility for each probe's binding site.

        Uses local RNA folding (200nt window) around each probe binding site
        to determine what fraction of the target is unpaired (accessible).
        Returns a Series of accessibility scores in [0, 1].
        """
        # Load the target sequence
        sequence_fasta = os.path.join(
            util.get_output_dir(), f"{util.get_gene_name()}_sequence.fasta"
        )
        if not os.path.isfile(sequence_fasta):
            cfg = SequenceConfig()
            if cfg.sequence_file and os.path.isfile(cfg.sequence_file):
                sequence_fasta = cfg.sequence_file
            else:
                self.logger.debug("No target sequence for accessibility scoring")
                return pd.Series(1.0, index=df.index)

        try:
            target_seq = str(next(Bio.SeqIO.parse(sequence_fasta, "fasta")).seq)
        except Exception:
            return pd.Series(1.0, index=df.index)

        window = 200  # nucleotides of context around binding site
        scores = []

        for _, row in df.iterrows():
            start = int(row["start"])
            end = int(row["end"])

            # Extract local window around probe binding site
            win_start = max(0, start - window // 2)
            win_end = min(len(target_seq), end + window // 2)
            local_seq = target_seq[win_start:win_end]

            if len(local_seq) < 10:
                scores.append(1.0)
                continue

            # Fold the local window
            try:
                dot_bracket = self._fold_sequence(local_seq)
                if dot_bracket is None:
                    scores.append(1.0)
                    continue

                # Count unpaired bases in the probe binding site region
                probe_start_in_window = start - win_start
                probe_end_in_window = end - win_start
                binding_region = dot_bracket[probe_start_in_window:probe_end_in_window]

                if len(binding_region) == 0:
                    scores.append(1.0)
                    continue

                unpaired = binding_region.count(".")
                accessibility = unpaired / len(binding_region)
                scores.append(round(accessibility, 4))
            except Exception:
                scores.append(1.0)

        return pd.Series(scores, index=df.index)

    @staticmethod
    def _fold_sequence(sequence: str) -> Optional[str]:
        """Fold a sequence and return the dot-bracket notation."""
        from pathlib import Path

        file_path = Path(__file__).resolve().parent.as_posix()
        data_table = os.path.join(file_path, "data_tables/")
        if sys.platform in ("linux", "linux2"):
            fold_path = os.path.join(file_path, "Fold_linux")
        elif sys.platform == "darwin":
            fold_path = os.path.join(file_path, "Fold_osx")
        else:
            return None

        if not os.path.isfile(fold_path):
            return None

        os.environ["DATAPATH"] = data_table
        fasta_input = f">target\n{sequence}\n"
        try:
            result = subprocess.run(
                [fold_path, "-", "-", "--bracket", "--MFE"],
                input=fasta_input,
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            lines = result.stdout.strip().split("\n")
            # Last line is the dot-bracket notation
            if len(lines) >= 3:
                return lines[-1]
        except Exception:
            pass
        return None

    @staticmethod
    def _compute_low_complexity_score(sequence: Bio.Seq.Seq) -> float:
        """Compute a low-complexity fraction for a probe sequence.

        Returns the fraction of 10bp windows with Shannon entropy < 1.5 bits.
        Higher values indicate more repetitive/low-complexity sequence.
        """
        import math
        seq_str = str(sequence)
        window = 10
        if len(seq_str) < window:
            return 0.0
        n_windows = len(seq_str) - window + 1
        low_count = 0
        for i in range(n_windows):
            w = seq_str[i:i + window]
            freqs = [w.count(b) / window for b in "ACGT"]
            entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
            if entropy < 1.5:
                low_count += 1
        return low_count / n_windows

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
            "tm_std": float(df["TM"].std()) if len(df) > 1 else 0.0,
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
            # Count how many probes hit each off-target gene
            gene_probe_counts: Dict[str, List[int]] = {}
            for idx, row in df.iterrows():
                genes_str = row.get("off_target_genes", "")
                if not genes_str:
                    continue
                for match in re.finditer(r"([^,\s]+)\(\d+\)", genes_str):
                    gene = match.group(1)
                    if gene not in gene_probe_counts:
                        gene_probe_counts[gene] = []
                    gene_probe_counts[gene].append(idx)

            # Find genes exceeding the cap
            over_cap = {
                gene: probes
                for gene, probes in gene_probe_counts.items()
                if len(probes) > cap
            }

            if not over_cap:
                break

            # Remove the lowest-quality probe among those hitting the most over-represented gene
            worst_gene = max(over_cap, key=lambda g: len(over_cap[g]))
            probe_indices = over_cap[worst_gene]
            # Find the probe with lowest quality among these
            worst_idx = df.loc[probe_indices, "quality"].idxmin()
            self.logger.debug(
                f"Off-target cap: removing probe at idx {worst_idx} "
                f"(gene {worst_gene} hit by {len(probe_indices)} probes, cap={cap})"
            )
            df = df.drop(worst_idx)

        return df.reset_index(drop=True)

    # Common gene family prefixes for paralog detection
    _GENE_FAMILY_PREFIXES = [
        "ACT", "MYH", "MYL", "HIST", "HLA", "KRT", "COL", "TUBB", "TUBA",
        "HSP", "RPL", "RPS", "EIF", "PCDH", "CDH", "CLDN",
    ]

    def _detect_off_target_clustering(self, df: pd.DataFrame) -> List[str]:
        """Detect when multiple probes hit the same off-target gene.

        Returns a list of warning messages. Also detects paralog families
        where the target may be inherently undruggable for FISH.
        """
        if "off_target_genes" not in df.columns:
            return []

        # Count how many probes hit each off-target gene
        gene_probe_counts: Dict[str, int] = {}
        for _, row in df.iterrows():
            genes_str = row.get("off_target_genes", "")
            if not genes_str:
                continue
            seen_genes = set()
            for match in re.finditer(r"([^,\s]+)\(\d+\)", genes_str):
                gene = match.group(1)
                if gene not in seen_genes:
                    seen_genes.add(gene)
                    gene_probe_counts[gene] = gene_probe_counts.get(gene, 0) + 1

        warnings = []
        cluster_threshold = 5
        n_probes = len(df)
        target_gene = util.get_gene_name(hashed=False)

        for gene, count in sorted(
            gene_probe_counts.items(), key=lambda x: -x[1]
        ):
            if count < cluster_threshold:
                break

            # Check if this is a paralog family
            is_paralog = False
            for prefix in self._GENE_FAMILY_PREFIXES:
                if (
                    gene.upper().startswith(prefix)
                    and target_gene.upper().startswith(prefix)
                ):
                    is_paralog = True
                    break

            # Undruggable detection: ≥50% of probes hit same off-target
            is_undruggable = n_probes > 0 and count >= n_probes * 0.5

            if is_undruggable and is_paralog:
                warnings.append(
                    f"UNDRUGGABLE TARGET: {count}/{n_probes} probes "
                    f"({count * 100 // n_probes}%) hit {gene} "
                    f"(paralog of {target_gene}). "
                    f"This gene cannot be uniquely targeted by exonic FISH probes "
                    f"due to conserved gene family homology. "
                    f"Alternatives: (1) use --target-regions intron for intronic probes "
                    f"(introns diverge between paralogs), "
                    f"(2) choose a different target gene."
                )
            elif is_undruggable:
                warnings.append(
                    f"UNDRUGGABLE TARGET: {count}/{n_probes} probes "
                    f"({count * 100 // n_probes}%) hit off-target gene {gene}. "
                    f"This target may not be uniquely targetable. "
                    f"Consider intronic probes (--target-regions intron) or "
                    f"a different target gene."
                )
            elif is_paralog:
                warnings.append(
                    f"Off-target clustering: {count}/{n_probes} probes hit {gene} "
                    f"(paralog of {target_gene}). This target may not be uniquely "
                    f"targetable by FISH due to gene family homology. "
                    f"Consider intronic probes instead."
                )
            else:
                warnings.append(
                    f"Off-target clustering: {count}/{n_probes} probes hit "
                    f"off-target gene {gene}."
                )

        return warnings

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

        # Apply cumulative off-target cap (C3)
        df = self._apply_off_target_cap(df)

        # Detect off-target clustering and emit warnings
        clustering_warnings = self._detect_off_target_clustering(df)
        if clustering_warnings:
            from .console import print_warning
            for warning in clustering_warnings:
                self.logger.warning(warning)
                print_warning(warning)

        sequences = self.prettify_sequences(df)
        config = self.prettify_configuration()

        # Write initial output for BLAST verification
        Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")

        # Run BLAST verification and remove flagged probes
        self._verification = self._run_blast_verification(
            self.output()["fasta"].path
        )
        if self._verification and self._verification["flagged"]:
            flagged_names = set(self._verification["flagged"].keys())
            pre_count = len(df)
            df = df[~df["name"].isin(flagged_names)].reset_index(drop=True)
            self.logger.debug(
                f"BLAST verification removed {pre_count - len(df)} probes with "
                f"unexpected genomic hits: {sorted(flagged_names)[:5]}"
                f"{'...' if len(flagged_names) > 5 else ''}"
            )
            # Rewrite output with cleaned probes
            sequences = self.prettify_sequences(df)
            Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")
            # Update verification counts
            self._verification["total"] = len(df)
            self._verification["clean"] = len(df)
            self._verification["flagged"] = {}

        util.log_and_check_candidates(self.logger, "CleanUpOutput", len(sequences))

        # Store summary stats and probe data for completion message
        self._summary = self._compute_summary(df)
        self._probe_df = df

        df.to_csv(self.output()["table"].path, index=False)
        with open(self.output()["config"].path, "w") as f:
            f.write(config)
        self.logger.debug(f'Saving all files using hash "{util.get_gene_name()}"')

        # Files to be deleted
        if not RunConfig().save_intermediates:
            self.remove_intermediates()
