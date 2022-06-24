"""Filter probes based on simple criteria."""

import logging
import multiprocessing
import os

import Bio.SeqIO
import Bio.SeqRecord
import Bio.SeqUtils
import Bio.SeqUtils.MeltingTemp
import luigi

from .config import GeneralConfig
from .config import ProbeConfig
from .generate_probes import GenerateAllProbes
from . import util


def get_melting_temp(
    sequence: str, na_concentration: float, formamide_concentration: float
) -> float:
    """Get melting temperature of candidate."""
    tm_raw = Bio.SeqUtils.MeltingTemp.Tm_NN(sequence, Na=na_concentration)
    melting_temp = Bio.SeqUtils.MeltingTemp.chem_correction(
        tm_raw, fmd=formamide_concentration
    )
    return melting_temp


def get_gc_content(sequence: str) -> float:
    """Get GC content of candidate."""
    return Bio.SeqUtils.GC(sequence)


def get_g_quadruplet_count(sequence: str) -> int:
    """Get number of G quadruplets in candidate."""
    return sequence.count("GGGG")


class BasicFiltering(luigi.Task):
    """Initial probe filtering based on melting temperature and GC content."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return GenerateAllProbes()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_basic.fasta")
        )

    def is_candidate_valid(self, candidate: Bio.SeqRecord.SeqRecord) -> bool:
        """Check if candidate matches basic criteria."""
        # TODO verify that min is less than max!
        sequence = candidate.seq
        gc_content = get_gc_content(sequence)
        g_quadruplets = get_g_quadruplet_count(sequence)
        melting_temp = get_melting_temp(
            sequence,
            ProbeConfig().na_concentration,
            ProbeConfig().formamide_concentration,
        )

        if (
            ProbeConfig().min_gc <= gc_content <= ProbeConfig().max_gc
            and ProbeConfig().min_tm <= melting_temp <= ProbeConfig().max_tm
            and g_quadruplets == 0
        ):
            return True
        return False

    def screen_sequences(self, sequences: list) -> list:
        """Filter candiates not matching criteria."""
        if len(sequences) > 1000:
            with multiprocessing.Pool(GeneralConfig().threads) as pool:
                valid = pool.map(self.is_candidate_valid, sequences)
        else:
            valid = [self.is_candidate_valid(candidate) for candidate in sequences]

        candidates = [seq for seq, valid in zip(sequences, valid) if valid]
        return candidates

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, "fasta"))
        candidates = self.screen_sequences(sequences)
        util.log_and_check_candidates(
            self.logger, "BasicFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, "fasta")
