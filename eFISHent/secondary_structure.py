"""Predict secondary structure of probes and filter them based on deltaG."""

from pathlib import Path
import logging
import multiprocessing
import os
import re
import subprocess
import sys
import tempfile

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .config import GeneralConfig
from .config import ProbeConfig
from .kmers import KMerFiltering


def get_free_energy(sequence: Bio.SeqRecord.SeqRecord) -> float:
    """Return the predicted free energy of the sequence.

    Using "Fold" as part of the RNAstructure package from the Mathews lab.

    Not ideal to include Fold_osx executable but easiest solution since
        there is no RNAstructure build on conda for MacOS and this is the only
        dependency not available.
    """
    file_path = Path(__file__).resolve().parent.as_posix()
    data_table = os.path.join(file_path, "data_tables/")
    if sys.platform == "linux" or sys.platform == "linux2":
        fold_path = "Fold"
    if sys.platform == "darwin":
        fold_path = os.path.join(file_path, "Fold_osx")

    os.environ["DATAPATH"] = data_table

    with tempfile.TemporaryDirectory() as tmp_dir:
        fname = os.path.join(tmp_dir, "input.fasta")
        Bio.SeqIO.write(sequence, fname, format="fasta")
        args_fold = [fold_path, fname, "-", "--bracket", "--MFE"]
        sec = subprocess.check_output(args_fold, stderr=subprocess.STDOUT).decode()

    # Format with secondary structure:
    # >ENERGY = -deltaG NAME DESCRIPTION\nSEQENCE\nDOT-BRACKET
    # Format without secondary structure:
    # >NAME DESCRIPTION\nSEQUENCE\n...
    search = re.search(r">ENERGY = (-\d+\.\d)\s\s([\w-]+)", sec)
    if search is not None:
        energy = float(search.groups()[0])
    else:
        energy = 0.0
    return energy


class SecondaryStructureFiltering(luigi.Task):
    """Predict free energy of secondary structures and filter them based on deltaG."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return KMerFiltering()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_sstruct.fasta")
        )

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, format="fasta"))

        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            free_energies = pool.map(get_free_energy, sequences)

        # Filter out candidates with low free energy
        candidates = [
            sequence
            for sequence, energy in zip(sequences, free_energies)
            if energy >= ProbeConfig().max_deltag
        ]

        util.log_and_check_candidates(
            self.logger, "SecondaryStructureFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
