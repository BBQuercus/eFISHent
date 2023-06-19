"""Filter probes based on simple criteria."""

from typing import List
import logging
import multiprocessing
import os

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.SeqUtils
import Bio.SeqUtils.MeltingTemp
import luigi

from . import util
from .config import GeneralConfig
from .config import ProbeConfig
from .generate_probes import GenerateAllProbes


def get_melting_temp(
    sequence: Bio.Seq.Seq, na_concentration: float, formamide_concentration: float
) -> float:
    """Get melting temperature of candidate assuming DNA/RNA hybrid."""
    rna_sequence = sequence.transcribe()
    tm_raw = Bio.SeqUtils.MeltingTemp.Tm_NN(
        rna_sequence,
        c_seq=rna_sequence.complement_rna(),
        nn_table=Bio.SeqUtils.MeltingTemp.R_DNA_NN1,
        Na=na_concentration,
    )
    melting_temp = Bio.SeqUtils.MeltingTemp.chem_correction(
        tm_raw, fmd=formamide_concentration
    )
    return float(melting_temp)


def get_gc_content(sequence: Bio.Seq.Seq) -> float:
    """Get GC content of candidate."""
    return float(Bio.SeqUtils.gc_fraction(sequence)) * 100


def get_g_quadruplet_count(sequence: Bio.Seq.Seq) -> int:
    """Get number of G quadruplets in candidate."""
    return int(sequence.count("GGGG"))


class BasicFiltering(luigi.Task):
    """Initial probe filtering based on melting temperature and GC content."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return GenerateAllProbes()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_basic.fasta")
        )

    def is_candidate_valid(
        self, candidate: Bio.SeqRecord.SeqRecord, config: luigi.Config
    ) -> bool:
        """Check if candidate matches basic criteria."""
        if config.min_tm > config.max_tm:
            raise ValueError("Max TM has to be greater or equal the min TM!")
        if config.min_gc > config.max_gc:
            raise ValueError("Max GC has to be greater or equal the min GC!")

        sequence = candidate.seq
        gc_content = get_gc_content(sequence)
        g_quadruplets = get_g_quadruplet_count(sequence)
        melting_temp = get_melting_temp(
            Bio.Seq.Seq(sequence),
            config.na_concentration,
            config.formamide_concentration,
        )

        if (
            config.min_gc <= gc_content <= config.max_gc
            and config.min_tm <= melting_temp <= config.max_tm
            and g_quadruplets == 0
        ):
            return True
        return False

    def _is_candidate_valid(self, candidate: Bio.SeqRecord.SeqRecord) -> bool:
        """Counteracting the really weird multiprocess class behavior."""
        return self.is_candidate_valid(candidate, self.config)

    def screen_sequences(
        self, sequences: List[Bio.SeqRecord.SeqRecord]
    ) -> List[Bio.SeqRecord.SeqRecord]:
        """Filter candiates not matching criteria."""
        if len(sequences) > 1000:
            self.config = ProbeConfig()
            with multiprocessing.Pool(GeneralConfig().threads) as pool:
                valid = pool.map(self._is_candidate_valid, sequences)
        else:
            valid = [
                self.is_candidate_valid(candidate, ProbeConfig())
                for candidate in sequences
            ]

        candidates = [seq for seq, valid in zip(sequences, valid) if valid]
        return candidates

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, "fasta"))
        candidates = self.screen_sequences(sequences)
        util.log_and_check_candidates(
            self.logger, "BasicFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
