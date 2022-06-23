import hashlib
import logging
import os
import re
import unicodedata

import luigi
import pandas as pd

from .config import GeneralConfig
from .config import SequenceConfig


FASTA_EXT = (".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa")


def get_output_dir() -> str:
    """Return the output directory."""
    output_dir = GeneralConfig().output_dir
    if output_dir != "None" and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        output_dir = os.getcwd()
    return output_dir


def hash_fn(key: bytearray) -> str:
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


def get_gene_name(hashed: bool = True) -> str:
    """Return the gene name without extension."""
    # Name based on gene fasta file
    sequence_file = SequenceConfig().sequence_file
    if sequence_file is not None and os.path.isfile(sequence_file):
        basename = secure_filename(os.path.splitext(os.path.basename(sequence_file))[0])
    # Using Ensembl ID
    elif SequenceConfig().ensemble_id:
        basename = secure_filename(SequenceConfig().ensemble_id)
    # Name based on NCBI parameters
    elif SequenceConfig().gene_name and SequenceConfig().organism_name:
        basename = secure_filename(
            f"{SequenceConfig().organism_name}_{SequenceConfig().gene_name}"
        )
    else:
        raise ValueError(
            "Could not determine gene name. "
            "Please specify the gene name and organism, or Ensembl ID, or pass in a sequence file."
        )

    if not hashed:
        return basename
    return f"{basename}_{create_config_hash(luigi.configuration.get_config())}"


def get_genome_name() -> str:
    """Return the genome name without extension."""
    if GeneralConfig().reference_genome is None:
        raise ValueError("Reference genome must be passed.")

    if not os.path.exists(
        GeneralConfig().reference_genome
    ) or not GeneralConfig().reference_genome.endswith(FASTA_EXT):
        raise ValueError("The passed reference genome must be a valid .fa file.")

    return os.path.abspath(os.path.splitext(GeneralConfig().reference_genome)[0])


def secure_filename(filename: str) -> str:
    """Pass a filename and return a secure version of it.
    This filename can then safely be stored on a regular file system and passed
    to `os.path.join`.  The filename returned is an ASCII only string for
    maximum portability.

    On windows system the function also makes sure that the file is not
    named after one of the special device files.

    >>> secure_filename("My cool movie.mov")
    'My_cool_movie.mov'
    >>> secure_filename("../../../etc/passwd")
    'etc_passwd'
    >>> secure_filename(u'i contain cool \xfcml\xe4uts.txt')
    'i_contain_cool_umlauts.txt'

    The function might return an empty filename.
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

    # on nt a couple of special files are present in each folder.  We
    # have to ensure that the target file is not such a filename.  In
    # this case we prepend an underline
    if (
        os.name == "nt"
        and filename
        and filename.split(".")[0].upper() in _windows_device_files
    ):
        filename = f"_{filename}"

    return filename


def log_and_check_candidates(
    logger: logging.Logger, name: str, count: int, count_prev: int = 0
) -> None:
    logger.info(f"Writing {count} (from {count_prev}) candidates in {name}.")

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


def create_data_table(sequences: list) -> pd.DataFrame:
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
