"""Create a list of candidate probes from a gene sequence."""

import logging
import os
from collections.abc import Generator, Iterator

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .config import ProbeConfig
from .prepare_sequence import PrepareSequence


def _compute_gc(seq: str) -> float:
    """Compute GC content as a fraction (0-1)."""
    if not seq:
        return 0.0
    gc = sum(1 for c in seq.upper() if c in "GC")
    return gc / len(seq)


def _preferred_length(gc_fraction: float, min_length: int, max_length: int) -> int:
    """Choose preferred probe length based on local GC content.

    High GC (>55%): prefer shorter probes to lower Tm
    Low GC (<45%): prefer longer probes to raise Tm
    Balanced (45-55%): use middle length
    """
    length_range = max_length - min_length
    if length_range < 2:
        return min_length

    if gc_fraction > 0.55:
        # High GC -> shorter probes
        return min_length
    elif gc_fraction < 0.45:
        # Low GC -> longer probes
        return max_length
    else:
        # Balanced GC -> middle length
        return min_length + length_range // 2


def _lengths_around_preferred(
    gc_fraction: float, min_length: int, max_length: int
) -> list[int]:
    """Get the GC-preferred probe length and its in-range neighbours, sorted."""
    preferred = _preferred_length(gc_fraction, min_length, max_length)
    lengths_to_try = {preferred}
    if preferred > min_length:
        lengths_to_try.add(preferred - 1)
    if preferred < max_length:
        lengths_to_try.add(preferred + 1)
    return sorted(lengths_to_try)


def _adaptive_windows(
    seq_str: str, min_length: int, max_length: int
) -> Iterator[tuple[int, int]]:
    """Yield (start, length) probe windows biased toward GC-optimal lengths.

    For each position, yields the preferred length ±1 to give the optimizer
    choices while biasing toward the GC-optimal length.
    """
    seq_length = len(seq_str)
    for start_pos in range(0, seq_length - min_length + 1):
        # Compute local GC in a window centered at this position
        window_start = max(0, start_pos - 10)
        window_end = min(seq_length, start_pos + max_length + 10)
        local_gc = _compute_gc(seq_str[window_start:window_end])

        for length in _lengths_around_preferred(local_gc, min_length, max_length):
            if start_pos + length <= seq_length:
                yield start_pos, length


def _all_windows(
    seq_length: int, min_length: int, max_length: int
) -> Iterator[tuple[int, int]]:
    """Yield (start, length) probe windows for all lengths at all positions."""
    for length in range(min_length, max_length + 1):
        for start_pos in range(0, seq_length - length + 1):
            yield start_pos, length


def create_candidate_probes_generator(
    sequence: Bio.SeqRecord.SeqRecord,
    min_length: int,
    max_length: int,
    adaptive: bool = False,
) -> Generator[Bio.SeqRecord.SeqRecord, None, None]:
    """Create a generator of all subsequences of sequence with the right lengths.

    This is a memory-efficient alternative to create_candidate_probes that yields
    probes one at a time instead of accumulating them all in memory.

    Args:
        sequence: The source sequence to generate probes from
        min_length: Minimum probe length
        max_length: Maximum probe length
        adaptive: If True, prefer lengths based on local GC content

    Yields:
        SeqRecord objects for each candidate probe
    """
    if min_length > max_length:
        raise ValueError(
            "Minimum probe length must be smaller or equal to maximum length. "
            f"{min_length} > {max_length}!"
        )
    if min_length >= len(sequence):
        raise ValueError(
            "Minimum probe length must be shorter than the sequence length. "
            f"{min_length} >= {len(sequence)}!"
        )

    if adaptive and max_length > min_length:
        windows = _adaptive_windows(str(sequence.seq).upper(), min_length, max_length)
    else:
        windows = _all_windows(len(sequence), min_length, max_length)

    for idx, (start_pos, length) in enumerate(windows, start=1):
        yield Bio.SeqRecord.SeqRecord(
            sequence.seq[start_pos : start_pos + length],
            id=f"candidate-{idx}-{start_pos}",
        )


class GenerateAllProbes(luigi.Task):
    """Create all possible probes in a gene given length ranges."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return PrepareSequence()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_all.fasta")
        )

    def create_candidate_probes(
        self, sequence: Bio.SeqRecord.SeqRecord, min_length: int, max_length: int
    ) -> list[Bio.SeqRecord.SeqRecord]:
        """Create a set of all subsequences of sequence with the right lengths."""
        if min_length > max_length:
            raise ValueError(
                "Minimum probe length must be smaller or equal to maximum length. "
                f"{min_length} > {max_length}!"
            )
        if min_length >= len(sequence):
            raise ValueError(
                "Minimum probe length must be shorter than the sequence length. "
                f"{min_length} >= {len(sequence)}!"
            )

        candidates = []
        idx = 1
        for length in range(min_length, max_length + 1):
            for start_pos in range(0, len(sequence) - length + 1):
                candidates.append(
                    Bio.SeqRecord.SeqRecord(
                        sequence.seq[start_pos : start_pos + length],
                        id=f"candidate-{idx}-{start_pos}",
                    )
                )
                idx += 1
        return candidates

    def run(self):
        util.log_stage_start(self.logger, "GenerateAllProbes")
        sequence = Bio.SeqIO.read(self.input().path, format="fasta")

        config = ProbeConfig()
        # Use generator for memory-efficient streaming to file
        count = 0
        with open(self.output().path, "w") as f:
            for candidate in create_candidate_probes_generator(
                sequence,
                config.min_length,
                config.max_length,
                adaptive=config.adaptive_length,
            ):
                Bio.SeqIO.write(candidate, f, format="fasta")
                count += 1

        util.log_and_check_candidates(self.logger, "GenerateAllProbes", count)
