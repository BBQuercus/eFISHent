import logging
import os
import re
import unicodedata

from config import GeneralConfig
from config import SequenceConfig


FASTA_EXT = (".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa")


def get_output_dir() -> str:
    """Return the output directory."""
    output_dir = GeneralConfig().output_dir
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        output_dir = os.getcwd()
    return output_dir


def get_gene_name() -> str:
    """Return the gene name without extension."""
    sequence_file = SequenceConfig().sequence_file
    if sequence_file is not None and os.path.isfile(sequence_file):
        return secure_filename(os.path.splitext(os.path.basename(sequence_file))[0])

    if SequenceConfig().gene_name and SequenceConfig().organism_name:
        return secure_filename(
            f"{SequenceConfig().organism_name}_{SequenceConfig().gene_name}"
        )

    raise ValueError(
        "Could not determine gene name. "
        "Please specify the gene name and organism or pass in a sequence file."
    )


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
    logger.debug(f"Writing {count} (from {count_prev}) candidates in {name}.")

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
