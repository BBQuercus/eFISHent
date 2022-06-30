"""Utility and getter functions."""

from typing import Any, List
import hashlib
import logging
import os
import re
import unicodedata

import Bio.SeqRecord
import luigi
import pandas as pd

from .config import GeneralConfig
from .config import SequenceConfig
from .constants import FASTA_EXT


def get_output_dir(config: luigi.Config = GeneralConfig) -> Any:
    """Return the output directory."""
    output_dir = config().output_dir
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not output_dir:
        output_dir = os.getcwd()
    return os.path.abspath(output_dir)


def hash_fn(key: str) -> str:
    """Deterministic hash function."""
    return hashlib.sha256(key.encode("utf-8")).hexdigest()[:10]


def create_config_hash(
    luigi_config: luigi.configuration.cfg_parser.LuigiConfigParser,
) -> str:
    """Create a hash from a luigi config parser."""
    return hash_fn(
        str(
            {
                section_name: dict(luigi_config[section_name])
                for section_name in luigi_config.sections()
            }
        )
    )


def get_gene_name(hashed: bool = True, config: luigi.Config = SequenceConfig) -> str:
    """Return the gene name without extension."""
    # Name based on gene fasta file
    sequence_file = config().sequence_file
    if sequence_file and os.path.isfile(sequence_file):
        basename = secure_filename(os.path.splitext(os.path.basename(sequence_file))[0])
    # Using Ensembl ID
    elif config().ensembl_id:
        basename = secure_filename(config().ensembl_id)
    # Name based on NCBI parameters
    elif config().gene_name and config().organism_name:
        basename = secure_filename(f"{config().organism_name}_{config().gene_name}")
    else:
        raise ValueError(
            "Could not determine gene name. "
            "Please specify the gene name and organism, or Ensembl ID, or pass in a sequence file."
        )

    if not hashed:
        return basename
    return f"{basename}_{create_config_hash(luigi.configuration.get_config())}"


def get_genome_name(config: luigi.Config = GeneralConfig) -> Any:
    """Return the genome name without extension."""
    if config().reference_genome is None:
        raise ValueError("Reference genome must be passed.")

    if not (
        os.path.exists(config().reference_genome)
        or os.path.splitext(config().reference_genome)[1] in FASTA_EXT
    ):
        raise ValueError("The passed reference genome must be a valid .fa file.")

    return os.path.abspath(os.path.splitext(config().reference_genome)[0])


def secure_filename(filename: str) -> str:
    """Pass a filename and return a secure version of it.

    This filename can then safely be stored on a regular file system and passed
    to `os.path.join`.  The filename returned is an ASCII only string for
    maximum portability.

    On windows system the function also makes sure that the file is not
    named after one of the special device files.
    """
    filename = unicodedata.normalize("NFKD", filename)
    filename = filename.encode("ascii", "ignore").decode("ascii")

    for sep in os.path.sep, os.path.altsep:
        if sep:
            filename = filename.replace(sep, " ")
    _filename_ascii_strip_re = re.compile(r"[^A-Za-z0-9_.-]")
    filename = str(_filename_ascii_strip_re.sub("", "_".join(filename.split()))).strip(
        "._"
    )
    return filename


def log_and_check_candidates(
    logger: logging.Logger, name: str, count: int, count_prev: int = 0
) -> None:
    """Log candidate count before/after filtering."""
    previous = f" (from {count_prev})" if count_prev else ""
    logger.info(f"Writing {count}{previous} candidates in {name}.")

    if count == 0:
        raise ValueError(
            f"No more probes remaining after {name}. "
            "Please check your gene sequence and the stringency of your parameters."
        )
    if count < 10:
        logger.warning(
            f"Only {count} candidates remain after {name}. "
            "Possibly check the stringency of your parameters."
        )


def create_data_table(sequences: List[Bio.SeqRecord.SeqRecord]) -> pd.DataFrame:
    """Create a data table from a list of sequences."""
    df = pd.DataFrame(
        list(map(lambda x: (x.id, len(x)), sequences)), columns=["name", "length"]
    )
    df.insert(1, "sequence", [str(seq.seq) for seq in sequences])
    df["start"] = (
        df["name"].str.extract(r"^candidate-\d+-(\d+)$", expand=True).astype(int)
    )
    df["end"] = df["start"] + df["length"]
    df = df.sort_values(by="start", ignore_index=True)
    return df


class UniCode:
    """Addition of fancy colors in help text."""

    ispos = os.name == "posix"
    blue = "\033[34m" if ispos else ""
    bold = "\033[1m" if ispos else ""
    cyan = "\033[36m" if ispos else ""
    green = "\033[32m" if ispos else ""
    magenta = "\033[35m" if ispos else ""
    red = "\033[31m" if ispos else ""
    yellow = "\033[33m" if ispos else ""
    end = "\033[0m" if ispos else ""
    dash = "\U0001F4A8"
    dna = "\U0001F9EC"
    error = "\U0001f6d1"
    fishing = "\U0001F3A3"
    happy = "\U0001F603"
    party = "\U0001F973"
    warn = "\u26A0"
