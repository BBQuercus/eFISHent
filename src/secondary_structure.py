"""
Predict secondary structure of probes and filter them based on deltaG.
"""

from typing import Tuple
import logging
import multiprocessing
import os
import re
import tempfile

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from config import GeneralConfig
from config import ProbeConfig
from kmers import KMerFiltering
import util


class SecondaryStructureFiltering(luigi.Task):
    """Predict free energy of secondary structures and filter them based on deltaG."""

    logger = logging.getLogger("luigi-interface")

    def requires(self):
        return KMerFiltering()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_sstruct.fasta")
        )

    @staticmethod
    def predict_secondary_structure(sequence: Bio.SeqRecord) -> str:
        with tempfile.TemporaryDirectory() as tmpdirname:
            fname = os.path.join(tmpdirname, "input.fasta")
            Bio.SeqIO.write(sequence, fname, "fasta")
            return os.popen(f"Fold {fname} - --bracket --MFE").read()

    def get_free_energy(self, sequence: Bio.SeqRecord) -> Tuple[str, float]:
        sec = self.predict_secondary_structure(sequence)
        # Format with secondary structure:
        # >ENERGY = -deltaG NAME DESCRIPTION\nSEQENCE\nDOT-BRACKET
        if search := re.search(r">ENERGY = (-\d+\.\d)\s\s([\w-]+)", sec):
            energy = float(search.groups()[0])
        # Format without secondary structure:
        # >NAME DESCRIPTION\nSEQUENCE\n...
        else:
            energy = 0.0
        return energy

    def run(self):
        sequences = list(Bio.SeqIO.parse(self.input().path, format="fasta"))

        # TODO inclue data tables in package
        os.environ["DATAPATH"] = "./data_tables/"
        with multiprocessing.Pool(GeneralConfig().threads) as pool:
            free_energies = pool.map(self.get_free_energy, sequences)
        self.logger.debug(f"Free energies: {free_energies}")

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
