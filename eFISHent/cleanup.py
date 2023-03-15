"""Clean up the output directory.

Remove the files that are not needed and prettify the kept output.
"""

from typing import List, Optional
import glob
import logging
import os
import sys

import luigi
import pandas as pd
import Bio.SeqIO

from . import util
from .alignment import AlignProbeCandidates
from .basic_filtering import get_gc_content
from .basic_filtering import get_melting_temp
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .kmers import BuildJellyfishIndex
from .kmers import get_max_kmer_count
from .optimization import OptimizeProbeCoverage
from .secondary_structure import get_free_energy


class CleanUpOutput(luigi.Task):
    """Clean up the output files and remove the intermediaries that are not needed."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_blast_required = (
            ProbeConfig().encode_count_table and SequenceConfig().is_endogenous
        )

    def requires(self):
        tasks = {
            "optimize": OptimizeProbeCoverage(),
            "jellyfish": BuildJellyfishIndex(),
        }
        if self.is_blast_required:
            tasks["alignment"] = AlignProbeCandidates()
        return tasks

    def output(self):
        return {
            name: luigi.LocalTarget(os.path.join(util.get_output_dir(), filename))
            for name, filename in [
                ("fasta", f"{util.get_gene_name()}.fasta"),
                ("table", f"{util.get_gene_name()}.csv"),
                ("config", f"{util.get_gene_name()}.txt"),
            ]
        }

    def prettify_table(
        self,
        sequences: List[Bio.SeqRecord.SeqRecord],
        basename: str,
        jellyfish_path: str,
        alignment_path: Optional[str] = None,
        config: luigi.Config = ProbeConfig,
    ) -> pd.DataFrame:
        """Create table with probe information."""
        # Create basic table with probe start/end positions and name
        df = util.create_data_table(sequences)
        sequences = sorted(sequences, key=lambda x: int(x.id.split("-")[-1]))

        # Add data columns
        df["GC"] = [round(get_gc_content(seq.seq), 2) for seq in sequences]
        df["TM"] = [
            round(
                get_melting_temp(
                    seq.seq, config().na_concentration, config().formamide_concentration
                ),
                2,
            )
            for seq in sequences
        ]
        df["deltaG"] = [get_free_energy(seq) for seq in sequences]
        df["kmers"] = [get_max_kmer_count(seq, jellyfish_path) for seq in sequences]
        if alignment_path is not None:
            df_counts = pd.read_csv(alignment_path)
            df_counts = df_counts.groupby("qname", as_index=False).size()
            df["count"] = pd.merge(
                df, df_counts, how="left", left_on="name", right_on="qname"
            )["size"].fillna(0)

        # Create new/clean names
        df["name"] = [f"{basename}-{idx + 1}" for idx in df.index]
        return df

    def prettify_sequences(self, df: pd.DataFrame) -> List[Bio.SeqRecord.SeqRecord]:
        """Clean up sequence id's and descriptions."""
        return [
            Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(seq), id=name, name=name, description=""
            )
            for seq, name in zip(df["sequence"], df["name"])
        ]

    def prettify_section(self, section: str) -> str:
        """Prettify a single configuration section."""
        pretty = [f"{section.replace('Config', '')} configuration:"]
        for key, value in self.config[section].items():
            if value and value != '""':
                pretty.append(f"  {key.replace('_', '-')}: {value}")
        return "\n".join(pretty)

    def prettify_configuration(self) -> str:
        """Create a pretty file with all set/default configuration."""
        self.config = luigi.configuration.get_config()
        pretty = [self.prettify_section(section) for section in self.config.sections()]
        command = " ".join(sys.argv[1:])
        pretty.append(f"Command:\n  efishent {command}")
        return "\n\n".join(pretty)

    def remove_intermediates(self) -> None:
        """Remove intermediary files no longer used."""
        intermediary_files = glob.glob(
            os.path.join(util.get_output_dir(), f"{util.get_gene_name()}_*")
        )
        self.logger.debug(f"Found file intermediates - {intermediary_files}")
        try:
            intermediary_files.remove(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_entrez.fasta"
                )
            )
        except ValueError:
            self.logger.debug("Didn't find entrez file...")
        for filename in intermediary_files:
            self.logger.debug(f"Removing {filename}")
            os.remove(filename)

    def run(self):
        sequences = list(
            Bio.SeqIO.parse(self.input()["optimize"]["probes"].path, "fasta")
        )
        df = self.prettify_table(
            sequences,
            basename=util.get_gene_name(hashed=False),
            jellyfish_path=self.input()["jellyfish"].path,
            alignment_path=self.input()["alignment"]["table"].path
            if self.is_blast_required
            else None,
        )
        sequences = self.prettify_sequences(df)
        config = self.prettify_configuration()

        util.log_and_check_candidates(self.logger, "CleanUpOutput", len(sequences))
        df.to_csv(self.output()["table"].path, index=False)
        Bio.SeqIO.write(sequences, self.output()["fasta"].path, format="fasta")
        with open(self.output()["config"].path, "w") as f:
            f.write(config)
        self.logger.info(f'Saving all files using hash "{util.get_gene_name()}"')

        # Files to be deleted
        if not RunConfig().save_intermediates:
            self.remove_intermediates()
