"""Index creating step for the alignment."""
import logging
import os
import subprocess
import warnings

import gtfparse
import luigi

from . import util
from .config import GeneralConfig


class BuildBowtieIndex(luigi.Task):
    """Create bowtie index for a given reference genome."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        return [
            luigi.LocalTarget(util.get_genome_name() + ext)
            for ext in [
                ".1.ebwt",
                ".2.ebwt",
                ".3.ebwt",
                ".4.ebwt",
                ".rev.1.ebwt",
                ".rev.2.ebwt",
            ]
        ]

    def run(self):
        args_bowtie = [
            "bowtie-build",
            os.path.abspath(GeneralConfig().reference_genome),
            util.get_genome_name(),
        ]
        self.logger.debug(f"Running bowtie with - {''.join(args_bowtie)}")
        subprocess.check_call(
            args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )
        self.logger.info("Finished building bowtie index.")


class PrepareAnnotationFile(luigi.Task):
    """Convert gtf annotation file to parquet.

    Doesn't change any data but greatly speeds up subsequent I/O (~10x increase).
    """

    logger = logging.getLogger("custom-logger")

    def output(self):
        return luigi.LocalTarget(GeneralConfig().reference_annotation + ".parq")

    def prepare_gtf_file(self, fname_input: str, fname_output: str) -> None:
        """Save GTF file as parquet."""
        warnings.filterwarnings("ignore")
        df = gtfparse.read_gtf(fname_input)
        warnings.filterwarnings("default")
        self.logger.debug("Read gtf file with gtfparse")
        df.to_parquet(fname_output)

    def run(self):
        self.prepare_gtf_file(GeneralConfig().reference_annotation, self.output().path)
        self.logger.info("Finished parsing GTF annotation file.")