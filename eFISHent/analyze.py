"""Analyze a probe set and save plots."""

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
import pysam

from . import util
from .alignment import AlignProbeCandidates
from .basic_filtering import get_g_quadruplet_count
from .basic_filtering import get_gc_content
from .basic_filtering import get_melting_temp
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .indexing import BuildBowtieIndex
from .kmers import BuildJellyfishIndex
from .kmers import get_max_kmer_count
from .prepare_sequence import PrepareSequence
from .secondary_structure import get_free_energy


class AnalyzeProbeset(luigi.Task):
    """Analysis of a probe set with summary stats."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        tasks = {
            "jellyfish": BuildJellyfishIndex(),
            "bowtie": BuildBowtieIndex(),
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
        return Bio.SeqIO.read(
            self.input()["sequence"].path, format="fasta"
        ).reverse_complement()

    @staticmethod
    def histplot(ax: plt.Axes, data: List[int], title: str, min_value: int):
        """Basic narrow count based histogram."""
        ax.set_title(title)
        bins = np.arange(min_value, max(data) + 4) + 0.5
        ax.hist(data, bins=bins)
        ax.set(xticks=bins[1:] - 0.5, ylabel="Count")

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
    def matrix(ax: plt.Axes, data: np.ndarray, title: str, fig):  # type: ignore
        """Image with colorbar in same axis."""
        ax.set_title(title)
        image = ax.imshow(data, vmin=0, vmax=1, cmap="coolwarm")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(image, cax=cax, orientation="vertical")

    def add_length(self, ax: plt.Axes):
        """Add probe lengths to axis."""
        lengths = list(map(len, self.sequences))
        self.histplot(ax, lengths, "Lengths", min(lengths) - 4)
        self.logger.debug("Added length to axis.")

    def add_melting_temperature(self, ax: plt.Axes):
        """Add probe melting temperatures to axis."""
        melting_temps = [
            get_melting_temp(
                seq.seq,
                ProbeConfig().na_concentration,
                ProbeConfig().formamide_concentration,
            )
            for seq in self.sequences
        ]
        self.boxplot(ax, melting_temps, "Melting temperatures", "ºC")
        self.logger.debug("Added melting temperature to axis.")

    def add_gc_content(self, ax: plt.Axes):
        """Add probe GC content to axis."""
        gc_content = [get_gc_content(seq.seq) for seq in self.sequences]
        self.boxplot(ax, gc_content, "GC Content", "%")
        self.logger.debug("Added GC content to axis.")

    def add_g_quadruplet(self, ax: plt.Axes):
        """Add number of G quadruplets in probes to axis."""
        quads = [get_g_quadruplet_count(seq.seq) for seq in self.sequences]
        self.histplot(ax, quads, "G quadruplet", -1)
        self.logger.debug("Added G quadruplets to axis.")

    def _get_max_kmer_count(self, seq):
        return get_max_kmer_count(seq, self.jellyfish_path)
        # return get_max_kmer_count(seq, "./tests/data/sacCer3_15.jf")

    def add_kmers(self, ax: plt.Axes):
        """Add highest number of sub-kmers in probes to axis."""
        self.jellyfish_path = self.input()["jellyfish"].path
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            kmers = pool.map(self._get_max_kmer_count, self.sequences)
        self.histplot(ax, kmers, "K-mer count (len)", -1)
        self.logger.debug("Added kmers to axis.")

    def add_free_energy(self, ax: plt.Axes):
        """Add predicted deltaG of probes to axis."""
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            deltag = pool.map(get_free_energy, self.sequences)
        self.boxplot(ax, deltag, "Free energy", "deltaG")
        self.logger.debug("Added free energy to axis.")

    def add_off_targets(self, ax: plt.Axes):
        """Add unique number of off targets to axis."""
        task = AlignProbeCandidates()
        task.fname_fasta = RunConfig().analyze_probeset
        task.fname_sam = os.path.splitext(task.fname_fasta)[0] + ".sam"
        task.fname_genome = util.get_genome_name()
        task.fname_gene = util.get_gene_name()
        task.is_endogenous = SequenceConfig().is_endogenous
        task.max_off_targets = 20

        task.align_probes(threads=GeneralConfig().threads)
        df_sam = task.parse_raw_pysam(pysam.view(task.fname_sam))  # type: ignore
        counts = df_sam.groupby("qname").size().values - 1
        self.histplot(ax, counts, "Off target count", -1)

        os.remove(task.fname_sam)
        os.remove(task.fname_fasta.rstrip("a") + "q")
        self.logger.debug("Added off targets to axis.")

    # TODO refactor to merge with optimization.is_binding
    def _get_similarity_score(self, seq1: Bio.Seq.Seq, seq2: Bio.Seq.Seq) -> float:
        alignments = [
            *Bio.pairwise2.align.globalxx(seq1, seq2.complement()),
            *Bio.pairwise2.align.globalxx(seq1, seq2.reverse_complement()),
        ]
        if not alignments:
            return 0.0
        max_score = max(map(lambda x: x.score, alignments))  # type: ignore
        return float(max_score / len(seq1))

    def add_binding_affinity(self, ax: plt.Axes, fig):
        """Add probe vs probe similarity scores (binding) to axis."""
        if len(self.sequences) >= 25:
            eta = len(self.sequences) ** 2 * 20 / 60_000
            self.logger.warning(
                f"Relatively large number of sequences ({len(self.sequences)}). "
                f"Determining sequence similarity will take ~{eta:.2f}mins."
            )
        matrix = np.array(
            [
                [self._get_similarity_score(x.seq, y.seq) for x in self.sequences]
                for y in self.sequences
            ]
        )  # type: ignore
        self.matrix(ax, matrix, "Binding affinity", fig)
        self.logger.debug("Added binding affinity to axis.")

    def add_probe_coverage(self, ax: plt.Axes):
        """Add coverage along gene to axis."""
        coverage = np.zeros(len(self.gene))

        for probe in self.sequences:
            idx = self.gene.seq.find(str(probe.seq))
            coverage[idx : idx + len(probe)] = 1

        ax.set_title("Gene coverage")
        ax.bar(range(len(self.gene)), coverage)
        ax.set(ylabel="Assigned", xlabel="Position")
        self.logger.debug("Added probe coverage to axis.")

    def build_figure(self):
        """Layout to build the figure."""
        fig = plt.figure(figsize=(15, 18))
        plt.suptitle(
            f"Analysis of probe set - {util.get_gene_name(False)}", y=1, fontsize=20
        )
        shape = (4, 3)

        # Row 1 - length, TM, GC
        ax = plt.subplot2grid(shape, (0, 0))
        self.add_length(ax)
        ax = plt.subplot2grid(shape, (0, 1))
        self.add_melting_temperature(ax)
        ax = plt.subplot2grid(shape, (0, 2))
        self.add_gc_content(ax)

        # Row 2 - G-quadruplet, kmer, deltaG
        ax = plt.subplot2grid(shape, (1, 0))
        self.add_g_quadruplet(ax)
        ax = plt.subplot2grid(shape, (1, 1))
        self.add_kmers(ax)
        ax = plt.subplot2grid(shape, (1, 2))
        self.add_free_energy(ax)

        # Row 3 - off target count, affinity
        ax = plt.subplot2grid(shape, (2, 0))
        self.add_off_targets(ax)
        ax = plt.subplot2grid(shape, (2, 1))
        self.add_binding_affinity(ax, fig)

        # Row 4 - coverage
        ax = plt.subplot2grid(shape, (3, 0), colspan=3)
        self.add_probe_coverage(ax)

        plt.tight_layout()
        plt.savefig(self.output().path)

    def run(self):
        self.sequences = list(
            Bio.SeqIO.parse(RunConfig().analyze_probeset, format="fasta")
        )
        self.build_figure()
        os.remove(self.input()["sequence"].path)
