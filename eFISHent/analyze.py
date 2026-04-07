"""Analyze a probe set and produce a detailed report with plots and metrics."""

from typing import List, Union
import logging
import multiprocessing
import os

from mpl_toolkits.axes_grid1 import make_axes_locatable
import Bio.Seq
import Bio.SeqIO
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from . import util
from .basic_filtering import compute_duplex_dg
from .basic_filtering import get_cpg_fraction
from .basic_filtering import get_gc_content
from .basic_filtering import get_g_quadruplet_count
from .basic_filtering import get_max_homopolymer_run
from .basic_filtering import get_melting_temp
from .basic_filtering import has_g_quadruplex
from .cleanup import CleanUpOutput
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .console import print_analysis_stage
from .indexing import BuildBowtieIndex
from .indexing import BuildBowtie2Index
from .kmers import BuildJellyfishIndex
from .kmers import get_max_kmer_counts_batch
from .optimization import is_binding
from .prepare_sequence import PrepareSequence
from .secondary_structure import get_free_energy


class AnalyzeProbeset(luigi.Task):
    """Analysis of a probe set with summary stats, plots, and CSV output."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        aligner = ProbeConfig().aligner
        if aligner == "bowtie2":
            index_task = BuildBowtie2Index()
        else:
            index_task = BuildBowtieIndex()
        tasks = {
            "jellyfish": BuildJellyfishIndex(),
            "bowtie": index_task,
            "sequence": PrepareSequence(),
        }
        return tasks

    def output(self):
        return luigi.LocalTarget(
            os.path.splitext(os.path.basename(RunConfig().analyze_probeset))[0]
            + "_analysis.pdf"
        )

    @property
    def gene(self):
        """Parse gene sequence file depending on parameters passed."""
        return Bio.SeqIO.read(self.input()["sequence"].path, format="fasta")

    # ── Plot helpers ──

    @staticmethod
    def histplot(ax: plt.Axes, data: List[int], title: str, min_value: int):
        """Basic narrow count based histogram."""
        ax.set_title(title)
        if not data or all(v == data[0] for v in data):
            ax.hist(data, bins=10)
        else:
            bins = np.arange(min_value, max(data) + 4) + 0.5
            ax.hist(data, bins=bins)
            ax.set(xticks=bins[1:] - 0.5)
        ax.set(ylabel="Count")

    @staticmethod
    def boxplot(ax: plt.Axes, data: List[Union[int, float]], title: str, ylabel: str):
        """Boxplot with kde-esque lines next to it."""
        ax.set_title(title)
        ax.boxplot(data)
        ax.eventplot(
            data,
            orientation="vertical",
            linewidths=0.2,
            lineoffsets=0.8,
            linelengths=0.1,
        )
        ax.set(xticks=[], ylabel=ylabel)

    @staticmethod
    def matrix(ax: plt.Axes, data: np.ndarray, title: str, fig):
        """Image with colorbar in same axis."""
        ax.set_title(title)
        image = ax.imshow(data, vmin=0, vmax=1, cmap="coolwarm")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(image, cax=cax, orientation="vertical")

    # ── Compute metrics ──

    def _compute_all_metrics(self) -> pd.DataFrame:
        """Compute all per-probe metrics and return as a DataFrame."""
        na = ProbeConfig().na_concentration
        formamide = ProbeConfig().formamide_concentration

        records = []
        for seq in self.sequences:
            s = seq.seq
            tm = get_melting_temp(s, na, formamide)
            gc = get_gc_content(s)
            cpg = get_cpg_fraction(s)
            dg = compute_duplex_dg(s, na, formamide)
            homopoly = get_max_homopolymer_run(s)
            g4 = has_g_quadruplex(s)
            g4_count = get_g_quadruplet_count(s)
            lc = CleanUpOutput._compute_low_complexity_score(s)

            records.append({
                "name": seq.id,
                "sequence": str(s),
                "length": len(s),
                "TM": round(tm, 2),
                "GC": round(gc, 2),
                "cpg_fraction": round(cpg, 4),
                "on_target_dg": round(dg, 2),
                "homopolymer_run": homopoly,
                "g_quadruplex": g4,
                "g_quadruplet_count": g4_count,
                "low_complexity": round(lc, 4),
            })

        df = pd.DataFrame(records)

        # K-mer counts
        jellyfish_path = self.input()["jellyfish"].path
        kmers = get_max_kmer_counts_batch(self.sequences, jellyfish_path)
        df["kmers"] = kmers

        # Secondary structure deltaG
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            deltag = pool.map(get_free_energy, self.sequences)
        df["deltaG"] = deltag

        # Probe positions on gene
        gene_seq = str(self.gene.seq)
        gene_rc = str(self.gene.seq.reverse_complement())
        starts, ends = [], []
        for seq in self.sequences:
            probe_seq = str(seq.seq)
            idx = gene_seq.find(probe_seq)
            if idx == -1:
                idx = gene_rc.find(probe_seq)
            if idx >= 0:
                starts.append(idx)
                ends.append(idx + len(seq))
            else:
                starts.append(-1)
                ends.append(-1)
        df["start"] = starts
        df["end"] = ends

        return df

    def _compute_off_targets(self, df: pd.DataFrame) -> pd.DataFrame:
        """Align probes and compute off-target counts."""
        from .alignment import AlignProbeCandidates

        task = AlignProbeCandidates()
        task.fname_fasta = RunConfig().analyze_probeset
        task.fname_sam = os.path.splitext(task.fname_fasta)[0] + "_analysis.sam"
        task.fname_genome = util.get_genome_name()
        task.is_endogenous = SequenceConfig().is_endogenous
        task.max_off_targets = 20

        aligner = ProbeConfig().aligner
        if aligner == "bowtie2":
            task.align_probes_bowtie2(threads=GeneralConfig().threads)
        else:
            task.align_probes(threads=GeneralConfig().threads)

        df_sam = task.parse_raw_pysam(pysam.view(task.fname_sam))
        counts = df_sam.groupby("qname").size()
        # Subtract 1 for the self-hit (endogenous)
        df["off_targets"] = df["name"].map(counts).fillna(0).astype(int) - 1
        df["off_targets"] = df["off_targets"].clip(lower=0)

        # Clean up temp files
        if os.path.isfile(task.fname_sam):
            os.remove(task.fname_sam)
        fastq = task.fname_fasta.rstrip("a") + "q"
        if os.path.isfile(fastq):
            os.remove(fastq)
        # Bowtie2 creates a clean FASTA copy
        bt2_fa = task.fname_fasta + ".bt2.fa"
        if os.path.isfile(bt2_fa):
            os.remove(bt2_fa)

        return df

    def _compute_binding_matrix(self) -> np.ndarray:
        """Compute probe-probe cross-hybridization matrix using k-mer overlap."""
        n = len(self.sequences)
        matrix = np.zeros((n, n))
        seqs = [str(s.seq) for s in self.sequences]

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    matrix[i, j] = 0.0
                    continue
                # Use k-mer overlap score (fraction of shared k-mers)
                seq1, seq2 = seqs[i], seqs[j]
                k = min(8, max(3, min(len(seq1), len(seq2)) - 1))
                kmers1 = {seq1[x:x + k] for x in range(len(seq1) - k + 1)}
                comp = str(Bio.Seq.Seq(seq2).complement())
                rc = str(Bio.Seq.Seq(seq2).reverse_complement())
                best_overlap = 0.0
                for target in (comp, rc):
                    kmers_t = {target[x:x + k] for x in range(len(target) - k + 1)}
                    denom = max(len(kmers1), len(kmers_t))
                    if denom > 0:
                        overlap = len(kmers1 & kmers_t) / denom
                        best_overlap = max(best_overlap, overlap)
                matrix[i, j] = best_overlap
                matrix[j, i] = best_overlap

        return matrix

    def _compute_quality_scores(self, df: pd.DataFrame) -> pd.Series:
        """Compute composite quality scores matching the main pipeline."""
        scores = pd.DataFrame(index=df.index)

        cfg = ProbeConfig()

        # Tm score: distance from median (narrower spread = better)
        if len(df) > 1:
            median_tm = df["TM"].median()
            tm_range = (cfg.max_tm - cfg.min_tm) / 2 if cfg.max_tm > cfg.min_tm else 5.0
            scores["tm"] = 1.0 - (df["TM"] - median_tm).abs().clip(upper=tm_range) / tm_range
        else:
            scores["tm"] = 1.0

        # GC score
        def _gc_score(gc: float) -> float:
            if 45.0 <= gc <= 52.0:
                return 1.0
            elif gc < 45.0:
                return max(0.0, (gc - 20.0) / 25.0)
            elif gc <= 55.0:
                return 1.0 - (gc - 52.0) / 10.0
            else:
                return max(0.0, 0.7 - (gc - 55.0) / 15.0)

        scores["gc"] = df["GC"].apply(_gc_score)

        # CpG score
        scores["cpg"] = (1.0 - df["cpg_fraction"].clip(upper=0.10) / 0.10).clip(lower=0.0)

        # deltaG score (secondary structure)
        dg_max = abs(cfg.max_deltag) if cfg.max_deltag < 0 else 10.0
        scores["dg"] = 1.0 - df["deltaG"].abs().clip(upper=dg_max) / dg_max

        # K-mer score
        kmer_max = max(cfg.max_kmers, 1)
        scores["kmer"] = 1.0 - df["kmers"].clip(upper=kmer_max) / kmer_max

        # Off-target score
        if "off_targets" in df.columns:
            ot_max = max(cfg.max_off_targets + 1, 1)
            scores["ot"] = 1.0 - df["off_targets"].clip(upper=ot_max) / ot_max
        else:
            scores["ot"] = 1.0

        # Low-complexity score
        scores["lc"] = 1.0 - df["low_complexity"].clip(upper=0.3) / 0.3

        # Binding ΔG score
        dg_floor = -60.0
        scores["binding"] = (df["on_target_dg"].clip(upper=0.0) / dg_floor).clip(
            lower=0.0, upper=1.0
        )

        # Clip all to [0, 1]
        for col in scores.columns:
            scores[col] = scores[col].clip(lower=0.0, upper=1.0)

        # Weighted composite
        quality = (
            scores["tm"] * 0.18
            + scores["gc"] * 0.13
            + scores["cpg"] * 0.05
            + scores["dg"] * 0.13
            + scores["ot"] * 0.23
            + scores["kmer"] * 0.09
            + scores["lc"] * 0.09
            + scores["binding"] * 0.10
        )
        return (quality * 100).round(1)

    # ── Summary output ──

    def _print_summary(self, df: pd.DataFrame):
        """Print a text summary to the terminal."""
        from .console import console

        n = len(df)
        gene_len = len(self.gene)
        mapped = df[df["start"] >= 0]

        console.print()
        console.print(f"  [bold]Probe Set Analysis[/bold] — {n} probes")
        console.print(f"  Gene length: {gene_len:,} nt")
        if len(mapped) > 0:
            coverage = mapped.apply(
                lambda r: r["end"] - r["start"], axis=1
            ).sum() / gene_len * 100
            console.print(f"  Gene coverage: {coverage:.1f}%")
            console.print(
                f"  Probes mapped to gene: {len(mapped)}/{n}"
            )
        console.print()

        # Per-metric summary table
        metrics = [
            ("Length", "length", "nt", 0),
            ("Tm", "TM", "°C", 1),
            ("GC", "GC", "%", 1),
            ("CpG fraction", "cpg_fraction", "", 3),
            ("ΔG (secondary)", "deltaG", "kcal/mol", 1),
            ("ΔG (binding)", "on_target_dg", "kcal/mol", 1),
            ("K-mer max", "kmers", "", 0),
            ("Homopolymer", "homopolymer_run", "nt", 0),
            ("Low complexity", "low_complexity", "", 3),
        ]

        console.print("  [dim]Metric              Median     Range[/dim]")
        for label, col, unit, dec in metrics:
            if col not in df.columns:
                continue
            vals = df[col]
            fmt = f".{dec}f"
            med = f"{vals.median():{fmt}}"
            lo = f"{vals.min():{fmt}}"
            hi = f"{vals.max():{fmt}}"
            suffix = f" {unit}" if unit else ""
            console.print(
                f"  {label:<20s} {med:>8s}   {lo} – {hi}{suffix}"
            )

        if "off_targets" in df.columns:
            ot = df["off_targets"]
            console.print(
                f"  {'Off-targets':<20s} {ot.median():>8.0f}   {ot.min():.0f} – {ot.max():.0f}"
            )

        if "quality" in df.columns:
            q = df["quality"]
            console.print(
                f"  {'Quality score':<20s} {q.median():>8.1f}   {q.min():.1f} – {q.max():.1f}"
            )

        # Flags
        flags = []
        g4_count = df["g_quadruplex"].sum()
        if g4_count > 0:
            flags.append(f"{g4_count} probe(s) contain G-quadruplex motifs")
        high_gc = (df["GC"] > 55).sum()
        if high_gc > 0:
            flags.append(f"{high_gc} probe(s) have GC > 55% (risk of non-specific binding)")
        high_lc = (df["low_complexity"] > 0.1).sum()
        if high_lc > 0:
            flags.append(f"{high_lc} probe(s) have elevated low-complexity score")
        if "off_targets" in df.columns:
            multi_ot = (df["off_targets"] > 0).sum()
            if multi_ot > 0:
                flags.append(f"{multi_ot} probe(s) have genome off-target hits")

        # Cross-hybridization
        binding_pairs = 0
        seqs = df["sequence"].values
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                if is_binding(seqs[i], seqs[j], 0.75):
                    binding_pairs += 1
        if binding_pairs > 0:
            flags.append(f"{binding_pairs} probe pair(s) may cross-hybridize (≥75% k-mer overlap)")

        if flags:
            console.print()
            console.print("  [bold yellow]Flags:[/bold yellow]")
            for f in flags:
                console.print(f"    [yellow]![/yellow] {f}")
        else:
            console.print()
            console.print("  [green]No issues detected[/green]")

        console.print()

    def _save_csv(self, df: pd.DataFrame):
        """Save per-probe metrics as CSV alongside the PDF."""
        csv_path = self.output().path.replace(".pdf", ".csv")
        # Select columns for output
        cols = [
            "name", "sequence", "length", "start", "end",
            "TM", "GC", "cpg_fraction", "on_target_dg",
            "deltaG", "kmers", "homopolymer_run",
            "g_quadruplex", "g_quadruplet_count", "low_complexity",
        ]
        if "off_targets" in df.columns:
            cols.append("off_targets")
        if "quality" in df.columns:
            cols.append("quality")
        out_cols = [c for c in cols if c in df.columns]
        df[out_cols].to_csv(csv_path, index=False)
        self.logger.info(f"Saved analysis CSV: {csv_path}")

    # ── Plot panels ──

    def _add_length(self, ax: plt.Axes, df: pd.DataFrame):
        lengths = df["length"].tolist()
        self.histplot(ax, lengths, "Lengths", min(lengths) - 4)

    def _add_melting_temperature(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["TM"].tolist(), "Melting temperatures", "°C")

    def _add_gc_content(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["GC"].tolist(), "GC Content", "%")

    def _add_cpg_fraction(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["cpg_fraction"].tolist(), "CpG Fraction", "")

    def _add_g_quadruplex(self, ax: plt.Axes, df: pd.DataFrame):
        quads = df["g_quadruplet_count"].tolist()
        self.histplot(ax, quads, "G-quadruplex count", -1)

    def _add_kmers(self, ax: plt.Axes, df: pd.DataFrame):
        self.histplot(ax, df["kmers"].tolist(), "K-mer count", -1)

    def _add_free_energy(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["deltaG"].tolist(), "Secondary structure ΔG", "kcal/mol")

    def _add_duplex_dg(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["on_target_dg"].tolist(), "On-target binding ΔG", "kcal/mol")

    def _add_off_targets(self, ax: plt.Axes, df: pd.DataFrame):
        if "off_targets" in df.columns:
            self.histplot(ax, df["off_targets"].tolist(), "Off-target count", -1)
        else:
            ax.set_title("Off-target count")
            ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)

    def _add_homopolymer(self, ax: plt.Axes, df: pd.DataFrame):
        self.histplot(ax, df["homopolymer_run"].tolist(), "Homopolymer run", 0)

    def _add_low_complexity(self, ax: plt.Axes, df: pd.DataFrame):
        self.boxplot(ax, df["low_complexity"].tolist(), "Low complexity score", "")

    def _add_binding_affinity(self, ax: plt.Axes, fig):
        if len(self.sequences) > 100:
            ax.set_title("Binding affinity")
            ax.text(
                0.5, 0.5, f"Skipped\n({len(self.sequences)} probes)",
                ha="center", va="center", transform=ax.transAxes,
            )
            return
        matrix = self._compute_binding_matrix()
        self.matrix(ax, matrix, "Binding affinity", fig)

    def _add_probe_coverage(self, ax: plt.Axes, df: pd.DataFrame):
        gene_len = len(self.gene)
        coverage = np.zeros(gene_len)
        mapped = df[df["start"] >= 0]

        for _, row in mapped.iterrows():
            s, e = int(row["start"]), int(row["end"])
            coverage[s:e] = 1

        ax.set_title("Gene coverage")
        ax.fill_between(range(gene_len), coverage, alpha=0.6, color="steelblue")

        # Mark probe positions
        for _, row in mapped.iterrows():
            mid = (row["start"] + row["end"]) / 2
            ax.axvline(mid, color="steelblue", alpha=0.3, linewidth=0.5)

        ax.set(ylabel="Coverage", xlabel="Position (nt)", ylim=(-0.1, 1.5))
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["No", "Yes"])

        # Add gap annotations
        if len(mapped) > 1:
            sorted_mapped = mapped.sort_values("start")
            prev_end = 0
            gaps = []
            for _, row in sorted_mapped.iterrows():
                gap = row["start"] - prev_end
                if gap > 100:
                    gaps.append((prev_end, int(row["start"]), gap))
                prev_end = max(prev_end, int(row["end"]))
            for gs, ge, size in gaps[:5]:  # Show max 5 largest gaps
                mid = (gs + ge) / 2
                ax.annotate(
                    f"{size}nt gap",
                    xy=(mid, 1.1), fontsize=6, ha="center", color="red", alpha=0.7,
                )

    def _add_quality_scores(self, ax: plt.Axes, df: pd.DataFrame):
        if "quality" not in df.columns:
            ax.set_title("Quality scores")
            ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)
            return
        self.boxplot(ax, df["quality"].tolist(), "Quality scores", "Score (0-100)")

    # ── Build figure ──

    def build_figure(self, df: pd.DataFrame):
        """Layout to build the figure."""
        fig = plt.figure(figsize=(18, 24))
        plt.suptitle(
            f"Probe Set Analysis — {util.get_gene_name(False)} ({len(self.sequences)} probes)",
            y=0.98, fontsize=16,
        )
        shape = (6, 3)

        # Row 1 — Length, Tm, GC
        self._add_length(plt.subplot2grid(shape, (0, 0)), df)
        self._add_melting_temperature(plt.subplot2grid(shape, (0, 1)), df)
        self._add_gc_content(plt.subplot2grid(shape, (0, 2)), df)

        # Row 2 — CpG, G-quadruplex, Homopolymer
        self._add_cpg_fraction(plt.subplot2grid(shape, (1, 0)), df)
        self._add_g_quadruplex(plt.subplot2grid(shape, (1, 1)), df)
        self._add_homopolymer(plt.subplot2grid(shape, (1, 2)), df)

        # Row 3 — K-mer, deltaG (secondary), deltaG (binding)
        self._add_kmers(plt.subplot2grid(shape, (2, 0)), df)
        self._add_free_energy(plt.subplot2grid(shape, (2, 1)), df)
        self._add_duplex_dg(plt.subplot2grid(shape, (2, 2)), df)

        # Row 4 — Off-targets, Binding affinity, Quality scores
        self._add_off_targets(plt.subplot2grid(shape, (3, 0)), df)
        self._add_binding_affinity(plt.subplot2grid(shape, (3, 1)), fig)
        self._add_quality_scores(plt.subplot2grid(shape, (3, 2)), df)

        # Row 5 — Low complexity (1 panel) + empty
        self._add_low_complexity(plt.subplot2grid(shape, (4, 0)), df)

        # Row 5-6 — Coverage (full width)
        ax_cov = plt.subplot2grid(shape, (5, 0), colspan=3)
        self._add_probe_coverage(ax_cov, df)

        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(self.output().path, dpi=150)
        plt.close()

    def run(self):
        total_steps = 7
        n_probes = len(list(Bio.SeqIO.parse(RunConfig().analyze_probeset, format="fasta")))
        print_analysis_stage(1, total_steps, f"Loading {n_probes} probes")
        self.sequences = list(
            Bio.SeqIO.parse(RunConfig().analyze_probeset, format="fasta")
        )

        # Compute all sequence-level metrics
        print_analysis_stage(2, total_steps, "Computing sequence metrics")
        df = self._compute_all_metrics()

        # Compute off-target counts
        print_analysis_stage(3, total_steps, "Aligning probes to genome")
        df = self._compute_off_targets(df)

        # Compute quality scores
        print_analysis_stage(4, total_steps, "Computing quality scores")
        df["quality"] = self._compute_quality_scores(df)

        # Print terminal summary
        print_analysis_stage(5, total_steps, "Summary")
        self._print_summary(df)

        # Save CSV
        print_analysis_stage(6, total_steps, "Saving CSV")
        self._save_csv(df)

        # Build and save PDF
        print_analysis_stage(7, total_steps, "Saving PDF")
        self.build_figure(df)

        os.remove(self.input()["sequence"].path)
