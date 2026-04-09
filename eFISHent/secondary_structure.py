"""Predict secondary structure of probes and filter them based on deltaG."""

from pathlib import Path
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys

import Bio.SeqIO
import Bio.SeqRecord
import luigi

from . import util
from .config import GeneralConfig
from .config import ProbeConfig
from .kmers import KMerFiltering

# Cached Fold binary path and whether to set DATAPATH
_fold_path: str = ""
_fold_datapath: str = ""


def _is_rnastructure_fold(path: str) -> bool:
    """Check whether a Fold binary is from the RNAstructure package."""
    try:
        out = subprocess.run(
            [path, "--version"],
            capture_output=True, text=True, timeout=5,
        )
        return "Mathews Lab" in (out.stdout + out.stderr)
    except Exception:
        return False


def _resolve_fold() -> None:
    """Find the best Fold binary and its DATAPATH (cached after first call)."""
    global _fold_path, _fold_datapath

    if _fold_path:
        return

    pkg_dir = Path(__file__).resolve().parent.as_posix()
    bundled_data = os.path.join(pkg_dir, "data_tables/")

    # Prefer system-installed Fold (from RNAstructure package) — it ships
    # with matching data tables, avoiding version mismatches with the
    # bundled binary.  Validate that the candidate is actually RNAstructure
    # (its --version output contains "Mathews Lab") to avoid picking up
    # unrelated binaries that happen to be named "Fold".
    system_fold = shutil.which("Fold")
    if system_fold and _is_rnastructure_fold(system_fold):
        _fold_path = system_fold
        # Keep existing DATAPATH if set (RNAstructure install sets it),
        # otherwise fall back to bundled data tables.
        _fold_datapath = os.environ.get("DATAPATH", bundled_data)
        return

    # Fall back to bundled binary
    if sys.platform in ("linux", "linux2"):
        _fold_path = os.path.join(pkg_dir, "Fold_linux")
    elif sys.platform == "darwin":
        _fold_path = os.path.join(pkg_dir, "Fold_osx")
    else:
        raise NotImplementedError(
            f"Platform '{sys.platform}' is not supported. "
            "eFISHent requires Linux or macOS for secondary structure prediction."
        )
    _fold_datapath = bundled_data


def get_free_energy(sequence: Bio.SeqRecord.SeqRecord) -> float:
    """Return the predicted free energy of the sequence.

    Using "Fold" as part of the RNAstructure package from the Mathews lab.

    Not ideal to include Fold_osx executable but easiest solution since
        there is no RNAstructure build for macOS and this is the only
        dependency not available.
    """
    _resolve_fold()
    os.environ["DATAPATH"] = _fold_datapath

    # Use stdin/stdout to avoid temp file creation overhead
    fasta_input = f">{sequence.id}\n{str(sequence.seq)}\n"
    args_fold = [_fold_path, "-", "-", "--bracket", "--MFE"]
    try:
        result = subprocess.run(
            args_fold,
            input=fasta_input,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        logger = logging.getLogger("custom-logger")
        logger.warning(
            f"Fold failed for {sequence.id} (exit code {e.returncode}): "
            f"{e.stderr.strip() or 'no error message'}. "
            "Assuming no secondary structure (deltaG = 0)."
        )
        return 0.0
    sec = result.stdout

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
        util.log_stage_start(self.logger, "SecondaryStructureFiltering")
        sequences = list(Bio.SeqIO.parse(self.input().path, format="fasta"))

        from .console import spinner

        with spinner(f"Predicting secondary structures ({len(sequences)} probes)..."):
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
