"""
Create a list of candidate probes from a gene sequence.
"""

import logging
import os

from tqdm import tqdm
import Bio.SeqIO
import Bio.SeqRecord
import Bio.SeqUtils
import Bio.SeqUtils.MeltingTemp
import luigi

from .config import ProbeConfig
from .prepare_sequence import PrepareSequence
from . import util


class GenerateAllProbes(luigi.Task):
    """Create all possible probes in a gene given length ranges."""

    logger = logging.getLogger("luigi-interface")

    def requires(self):
        return PrepareSequence()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_all.fasta")
        )

    def run(self):
        sequence = Bio.SeqIO.read(self.input().path, format="fasta")

        candidates = []
        idx = 1
        for length in tqdm(
            range(ProbeConfig().min_length, ProbeConfig().max_length + 1)
        ):
            for i in range(0, len(sequence.seq) - length + 1):
                candidates.append(
                    Bio.SeqRecord.SeqRecord(
                        sequence.seq[i : i + length], id=f"candidate-{idx}-{i}"
                    )
                )
                idx += 1

        util.log_and_check_candidates(self.logger, "GenerateAllProbes", len(candidates))
        Bio.SeqIO.write(candidates, self.output().path, "fasta")
