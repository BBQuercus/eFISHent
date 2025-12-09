"""Create jellyfish indices of the reference genome.

Filter probes based on jellyfish indices.
"""

import logging
import os
import subprocess
import tempfile
from typing import List

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .alignment import AlignProbeCandidates
from .config import GeneralConfig
from .config import ProbeConfig


def get_max_kmer_count(sequence: Bio.SeqRecord.SeqRecord, jellyfish_path: str) -> int:
    """Count kmers in a sequence and return the maximum count.

    Uses file-based input to jellyfish for handling of all k-mers.

    Args:
        sequence: SeqRecord object to analyze
        jellyfish_path: Path to the jellyfish index file

    Returns:
        Maximum kmer count found in the sequence
    """
    kmer_length = ProbeConfig().kmer_length
    if len(sequence.seq) < kmer_length:
        raise ValueError(
            f"Probe length must be larger than kmer length ({kmer_length})."
        )

    sub_kmers = [
        str(sequence.seq[i : i + kmer_length])
        for i in range(len(sequence) - kmer_length + 1)
    ]

    # Write k-mers to temp file in FASTA format (required by jellyfish -s)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        for i, kmer in enumerate(sub_kmers):
            f.write(f">kmer{i}\n{kmer}\n")
        temp_path = f.name

    try:
        args_jellyfish = ["jellyfish", "query", jellyfish_path, "-s", temp_path]
        kmer_counts_raw = subprocess.check_output(
            args_jellyfish, stderr=subprocess.STDOUT
        ).decode()
    finally:
        os.unlink(temp_path)

    # Format of jellyfish output: "KMER COUNT\n"
    kmer_counts = [
        int(line.split(" ")[1]) for line in kmer_counts_raw.strip().split("\n") if line
    ]
    return max(kmer_counts) if kmer_counts else 0


def get_max_kmer_counts_batch(
    sequences: List[Bio.SeqRecord.SeqRecord], jellyfish_path: str
) -> List[int]:
    """Count kmers for multiple sequences in a single jellyfish query.

    This is more efficient than calling get_max_kmer_count for each sequence
    individually, as it reduces subprocess overhead by batching all k-mers
    into a single jellyfish query using a temp file.

    Args:
        sequences: List of SeqRecord objects to analyze
        jellyfish_path: Path to the jellyfish index file

    Returns:
        List of max kmer counts, one per input sequence
    """
    if not sequences:
        return []

    kmer_length = ProbeConfig().kmer_length

    # Validate all sequences have sufficient length
    for seq in sequences:
        if len(seq.seq) < kmer_length:
            raise ValueError(
                f"Probe length must be larger than kmer length ({kmer_length})."
            )

    # Extract all kmers from all sequences, tracking which sequence each belongs to
    all_kmers = []
    kmer_counts_per_seq = []  # Number of kmers per sequence

    for sequence in sequences:
        seq_kmers = [
            str(sequence.seq[i : i + kmer_length])
            for i in range(len(sequence) - kmer_length + 1)
        ]
        all_kmers.extend(seq_kmers)
        kmer_counts_per_seq.append(len(seq_kmers))

    # Write all kmers to temp file in FASTA format (required by jellyfish -s)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        for i, kmer in enumerate(all_kmers):
            f.write(f">kmer{i}\n{kmer}\n")
        temp_path = f.name

    try:
        args_jellyfish = ["jellyfish", "query", jellyfish_path, "-s", temp_path]
        kmer_counts_raw = subprocess.check_output(
            args_jellyfish, stderr=subprocess.STDOUT
        ).decode()
    finally:
        os.unlink(temp_path)

    # Parse all counts (format: "KMER COUNT\n")
    all_counts = [
        int(line.split(" ")[1]) for line in kmer_counts_raw.strip().split("\n") if line
    ]

    # Split counts back per sequence and find max for each
    result = []
    idx = 0
    for count in kmer_counts_per_seq:
        seq_counts = all_counts[idx : idx + count]
        result.append(max(seq_counts) if seq_counts else 0)
        idx += count

    return result


class BuildJellyfishIndex(luigi.Task):
    """Index building task for jellyfish kmers."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(
            f"{util.get_genome_name()}_{ProbeConfig().kmer_length}.jf"
        )

    def run(self):
        util.log_stage_start(self.logger, "BuildJellyfishIndex")
        if ProbeConfig().max_kmers <= 2:
            self.logger.warning(
                f"{util.UniCode.warn} Jellyfish index will be created but not used because max_kmers <= 2."
            )

        args_jellyfish = [
            "jellyfish",
            "count",
            "--mer-len",
            str(ProbeConfig().kmer_length),
            "--out-counter-len",
            "1",
            "--lower-count",
            "2",
            "--size",
            "100M",
            "--threads",
            str(GeneralConfig().threads),
            "--output",
            self.output().path,
            GeneralConfig().reference_genome,
        ]
        self.logger.debug(f"Running jellyfish with - {' '.join(args_jellyfish)}")
        from .console import spinner

        with spinner("Building k-mer index..."):
            subprocess.check_call(args_jellyfish)
        self.logger.info(f"Finished building kmer ({ProbeConfig().kmer_length}) index.")


class KMerFiltering(luigi.Task):
    """Filter probes containing too many common short kmers."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return {"jellyfish": BuildJellyfishIndex(), "probes": AlignProbeCandidates()}

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_kmer.fasta")
        )

    def run(self):
        util.log_stage_start(self.logger, "KMerFiltering")
        sequences = list(
            Bio.SeqIO.parse(self.input()["probes"]["fasta"].path, format="fasta")
        )

        jellyfish_path = self.input()["jellyfish"].path
        counts = get_max_kmer_counts_batch(sequences, jellyfish_path)

        candidates = [
            seq
            for seq, count in zip(sequences, counts)
            if count <= ProbeConfig().max_kmers
        ]

        util.log_and_check_candidates(
            self.logger, "KMerFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
