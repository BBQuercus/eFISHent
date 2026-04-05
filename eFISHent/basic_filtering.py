"""Filter probes based on simple criteria."""

from typing import List
import logging
import math
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


def get_max_homopolymer_run(sequence: Bio.Seq.Seq) -> int:
    """Get the length of the longest homopolymer run in any base.

    Approach matches OligoMiner's blockParse.py prohibited sequences
    (Beliveau et al., 2018, MIT license).
    """
    seq_str = str(sequence)
    if not seq_str:
        return 0
    max_run = 1
    current_run = 1
    for i in range(1, len(seq_str)):
        if seq_str[i] == seq_str[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def get_dinucleotide_repeat_count(sequence: Bio.Seq.Seq, min_repeats: int = 4) -> int:
    """Count distinct dinucleotide repeat motifs present (e.g., ATATAT)."""
    seq_str = str(sequence)
    count = 0
    seen = set()
    for i in range(len(seq_str) - 1):
        di = seq_str[i : i + 2]
        if di not in seen and di[0] != di[1]:
            seen.add(di)
            if (di * min_repeats) in seq_str:
                count += 1
    return count


def get_cpg_fraction(sequence: Bio.Seq.Seq) -> float:
    """Get the fraction of CpG dinucleotides in the sequence.

    CpG frequency is significantly elevated in probes causing nuclear/nucleolar
    background (p=0.004, n=84 probe sets), likely from cross-hybridization with
    GC-rich rRNA (~60-67% GC).
    """
    seq_str = str(sequence).upper()
    if len(seq_str) < 2:
        return 0.0
    cpg_count = sum(1 for i in range(len(seq_str) - 1) if seq_str[i:i+2] == "CG")
    return cpg_count / (len(seq_str) - 1)


def has_low_complexity(
    sequence: Bio.Seq.Seq, window: int = 10, threshold: float = 1.0
) -> bool:
    """Check if any window of the sequence has Shannon entropy < threshold bits."""
    seq_str = str(sequence)
    for i in range(len(seq_str) - window + 1):
        w = seq_str[i : i + window]
        freqs = [w.count(b) / window for b in "ACGT"]
        entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
        if entropy < threshold:
            return True
    return False


def compute_off_target_tm(
    probe_seq: str,
    target_seq: str,
    na_concentration: float,
    formamide_concentration: float,
) -> float:
    """Estimate Tm for a probe binding to an off-target genomic sequence.

    Computes perfect-match RNA/DNA hybrid Tm using nearest-neighbor model,
    then subtracts ~2.5°C per mismatch (empirical penalty from SantaLucia 1998).
    Returns 0.0 on any error (e.g. non-standard bases).
    """
    min_len = min(len(probe_seq), len(target_seq))
    probe_trimmed = probe_seq[:min_len]
    target_trimmed = target_seq[:min_len]
    mismatches = sum(
        1 for a, b in zip(probe_trimmed, target_trimmed) if a != b
    )
    try:
        rna_probe = Bio.Seq.Seq(probe_trimmed).transcribe()
        perfect_complement = rna_probe.complement_rna()
        tm_raw = Bio.SeqUtils.MeltingTemp.Tm_NN(
            rna_probe,
            c_seq=perfect_complement,
            nn_table=Bio.SeqUtils.MeltingTemp.R_DNA_NN1,
            Na=na_concentration,
        )
        tm = Bio.SeqUtils.MeltingTemp.chem_correction(
            tm_raw, fmd=formamide_concentration
        )
        tm -= mismatches * 2.5
        return float(tm)
    except Exception:
        return 0.0


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
        melting_temp = get_melting_temp(
            Bio.Seq.Seq(sequence),
            config.na_concentration,
            config.formamide_concentration,
        )

        if not (config.min_gc <= gc_content <= config.max_gc):
            return False
        if not (config.min_tm <= melting_temp <= config.max_tm):
            return False

        # Homopolymer filtering
        max_homopolymer = getattr(config, "max_homopolymer_length", 0)
        if max_homopolymer > 0:
            if get_max_homopolymer_run(sequence) >= max_homopolymer:
                return False
        else:
            # Legacy fallback: only filter GGGG
            if get_g_quadruplet_count(sequence) > 0:
                return False

        # CpG depletion filter (optional)
        max_cpg = getattr(config, "max_cpg_fraction", 0.0)
        if max_cpg > 0:
            if get_cpg_fraction(sequence) > max_cpg:
                return False

        # Low-complexity filter (optional)
        if getattr(config, "filter_low_complexity", False):
            if has_low_complexity(sequence) or get_dinucleotide_repeat_count(sequence) > 0:
                return False

        return True

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
        util.log_stage_start(self.logger, "BasicFiltering")
        sequences = list(Bio.SeqIO.parse(self.input().path, "fasta"))
        candidates = self.screen_sequences(sequences)
        util.log_and_check_candidates(
            self.logger, "BasicFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
