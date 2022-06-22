"""
Create jellyfish indices of the reference genome.
Filter probes based on jellyfish indices.
"""

from typing import Optional
import logging
import multiprocessing
import os

import luigi
import Bio.SeqIO
import Bio.SeqRecord

from config import GeneralConfig
from config import ProbeConfig
from alignment import AlignProbeCandidates
import util


def get_max_kmer_count(sequence: Bio.SeqRecord, jellyfish_path: str) -> int:
    """Count kmers in a sequence.

    Only keep the sequence if it has less than max_kmers kmers.
    """

    sub_kmers = [
        str(sequence.seq[i : i + ProbeConfig().kmer_length])
        for i in range(len(sequence) - ProbeConfig().kmer_length + 1)
    ]
    kmer_counts_raw = os.popen(
        f"jellyfish query\
            {jellyfish_path}\
            {' '.join(sub_kmers)}"
    ).read()

    # Format of jellyfish output: "kmer1 count1\nkmer2 count2\n..."
    kmer_counts = list(
        map(lambda x: int(x.split(" ")[1]), kmer_counts_raw.split("\n")[:-1])
    )
    return max(kmer_counts)


class BuildJellyfishIndex(luigi.Task):
    """Index building task for jellyfish kmers."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(
            f"{util.get_genome_name()}_{ProbeConfig().kmer_length}.jf"
        )

    def run(self):
        if ProbeConfig().max_kmers <= 2:
            self.logger.warning(
                f"{util.UniCode.warn} Jellyfish index will be created but not used because max_kmers <= 2."
            )

        os.system(
            f"jellyfish count\
                --mer-len {ProbeConfig().kmer_length}\
                --out-counter-len 1\
                --lower-count 2\
                --size 100M\
                --threads {GeneralConfig().threads}\
                --output {self.output().path}\
                {GeneralConfig().reference_genome}"
        )


class KMerFiltering(luigi.Task):
    """Filter probes containing too many common short kmers."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return {"jellyfish": BuildJellyfishIndex(), "probes": AlignProbeCandidates()}

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_kmer.fasta")
        )

    def _get_kmer_count(self, sequence: Bio.SeqRecord) -> int:
        return get_max_kmer_count(sequence, self.input()["jellyfish"].path)

    def run(self):
        sequences = list(
            Bio.SeqIO.parse(self.input()["probes"]["fasta"].path, format="fasta")
        )

        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            counts = pool.map(self._get_kmer_count, sequences)

        candidates = [
            seq
            for seq, count in zip(sequences, counts)
            if count <= ProbeConfig().max_kmers
        ]

        util.log_and_check_candidates(
            self.logger, "KMerFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
