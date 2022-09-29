"""Create a list of candidate probes from a gene sequence."""

import logging
import os
from typing import List

import Bio.SeqIO
import Bio.SeqRecord
import Bio.SeqUtils
import Bio.SeqUtils.MeltingTemp
import luigi

from . import util
from .config import ProbeConfig
from .prepare_sequence import PrepareSequence


class GenerateAllProbes(luigi.Task):
    """Create all possible probes in a gene given length ranges."""

    logger = logging.getLogger("luigi-interface")

    def requires(self):
        return PrepareSequence()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_all.fasta")
        )

    def create_candidate_probes(
        self, sequence: Bio.SeqRecord.SeqRecord, min_length: int, max_length: int
    ) -> List[Bio.SeqRecord.SeqRecord]:
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
        sequence = Bio.SeqIO.read(self.input().path, format="fasta")
        candidates = self.create_candidate_probes(
            sequence, ProbeConfig().min_length, ProbeConfig().max_length
        )
        util.log_and_check_candidates(self.logger, "GenerateAllProbes", len(candidates))
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
