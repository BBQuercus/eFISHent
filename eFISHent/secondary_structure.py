"""
Predict secondary structure of probes and filter them based on deltaG.
"""

import logging
import multiprocessing
import os
import re
import subprocess
import tempfile

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from .config import GeneralConfig
from .config import ProbeConfig
from .kmers import KMerFiltering
from . import util


def get_free_energy(sequence: Bio.SeqRecord) -> float:
    """Return the predicted free energy of the sequence."""
    os.environ["DATAPATH"] = "./data_tables/"

    with tempfile.TemporaryDirectory() as tmpdirname:
        fname = os.path.join(tmpdirname, "input.fasta")
        Bio.SeqIO.write(sequence, fname, "fasta")
        args_fold = ["Fold", fname, "-", "--bracket", "--MFE"]
        sec = subprocess.check_output(args_fold, stderr=subprocess.STDOUT).decode()

    # Format with secondary structure:
    # >ENERGY = -deltaG NAME DESCRIPTION\nSEQENCE\nDOT-BRACKET
    if search := re.search(r">ENERGY = (-\d+\.\d)\s\s([\w-]+)", sec):
        energy = float(search.groups()[0])
    # Format without secondary structure:
    # >NAME DESCRIPTION\nSEQUENCE\n...
    else:
        energy = 0.0
    return energy


class SecondaryStructureFiltering(luigi.Task):
    """Predict free energy of secondary structures and filter them based on deltaG."""

    logger = logging.getLogger("luigi-interface")

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
            if energy >= ProbeConfig().max_deltaG
        ]

        util.log_and_check_candidates(
            self.logger, "SecondaryStructureFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
