"""Filter probes against ribosomal DNA and other high-copy dangerous sequences.

The 45S rDNA repeat unit (U13369.1, 42,999 bp) is absent from standard genome
assemblies (GRCh38, GRCh37) because the ~300-400 tandem copies reside on
unassembled acrocentric chromosome short arms. A single probe with even 2-3
mismatches to 45S rDNA produces overwhelming nucleolar signal in FISH experiments.

This filter is mandatory and runs regardless of endogenous/exogenous status,
because neither Bowtie2 genome alignment nor GTF-based rRNA annotation can
detect these off-targets.
"""

import logging
import os
import shutil
import subprocess
import tempfile

import Bio.SeqIO
import luigi

from . import util
from .basic_filtering import BasicFiltering
from .config import ProbeConfig


# Bundled reference sequences shipped with the package
_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RDNA_45S_FASTA = os.path.join(_DATA_DIR, "rDNA_45S_U13369.fasta")
ALPHA_SAT_FASTA = os.path.join(_DATA_DIR, "alpha_satellite_consensus.fasta")
IMAGING_RISK_FASTA = os.path.join(_DATA_DIR, "imaging_risk_rnas.fasta")


class FilterRibosomalRNA(luigi.Task):
    """BLAST probes against dangerous high-copy sequences.

    Screens against:
    - 45S rDNA (U13369.1): ~300-400 copies, absent from GRCh38, produces
      nucleolar signal
    - Alpha satellite consensus: ~millions of copies at centromeres
    - Imaging-risk RNA panel: mitochondrial rRNAs, U1/U2/U6 snRNAs,
      7SL/SRP RNA, 7SK RNA, and representative tRNAs — abundant RNAs
      that create concentrated background when hit by even a few probes

    This filter is mandatory and independent of user-provided references.
    """

    logger = logging.getLogger("custom-logger")

    def requires(self):
        return BasicFiltering()

    def output(self):
        return luigi.LocalTarget(
            os.path.join(
                util.get_output_dir(),
                f"{util.get_gene_name()}_rdna_filtered.fasta",
            )
        )

    def _build_blast_db(self, fasta_path: str, db_path: str) -> bool:
        """Build a temporary BLAST nucleotide database."""
        if not shutil.which("makeblastdb"):
            return False
        try:
            subprocess.check_call(
                ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", db_path],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
            return True
        except subprocess.CalledProcessError:
            self.logger.warning(f"Failed to build BLAST DB from {fasta_path}")
            return False

    def _blast_probes(
        self, probe_fasta: str, db_path: str, blast_out: str, threads: int
    ) -> None:
        """BLAST probes against a small reference database."""
        subprocess.check_call(
            [
                "blastn",
                "-task", "blastn",
                "-query", probe_fasta,
                "-db", db_path,
                "-evalue", "1000",
                "-word_size", "7",
                "-gapopen", "5",
                "-gapextend", "2",
                "-reward", "1",
                "-penalty", "-3",
                "-dust", "no",
                "-num_alignments", "100",
                "-num_threads", str(threads),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                           "qstart qend sstart send evalue bitscore qlen",
                "-out", blast_out,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

    def _find_dangerous_probes(self, blast_out: str) -> set:
        """Parse BLAST output and find probes with dangerous near-matches.

        Rejection criteria: any probe with <=3 mismatches across >=85% of
        its length against any position in the reference. This covers the
        entire pre-rRNA precursor including ETS and ITS regions.
        """
        import pandas as pd

        if not os.path.isfile(blast_out) or os.path.getsize(blast_out) == 0:
            return set()

        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen",
        ]
        df = pd.read_csv(blast_out, sep="\t", header=None, names=columns)

        dangerous = set()
        for _, row in df.iterrows():
            probe_len = row["qlen"]
            effective_len = row["length"] - row["gapopen"]
            mismatches = row["mismatch"]

            # Reject if alignment covers >=85% of probe and has <=3 mismatches
            if effective_len >= 0.85 * probe_len and mismatches <= 3:
                dangerous.add(row["qseqid"])

        return dangerous

    def run(self):
        from .config import GeneralConfig
        from .console import spinner

        util.log_stage_start(self.logger, "FilterRibosomalRNA")
        probe_fasta = self.input().path
        sequences = list(Bio.SeqIO.parse(probe_fasta, "fasta"))

        if not shutil.which("blastn") or not shutil.which("makeblastdb"):
            self.logger.warning(
                "BLAST+ not found — skipping rDNA/satellite filter. "
                "Install BLAST+ to enable this safety check."
            )
            Bio.SeqIO.write(sequences, self.output().path, format="fasta")
            return

        threads = GeneralConfig().threads
        all_dangerous = set()

        # Collect reference FASTAs to screen against
        custom_rdna = ProbeConfig().custom_rdna_fasta
        references = []

        if custom_rdna and os.path.isfile(custom_rdna):
            references.append(("custom rDNA", custom_rdna))
        elif os.path.isfile(RDNA_45S_FASTA):
            references.append(("45S rDNA (U13369.1)", RDNA_45S_FASTA))

        if os.path.isfile(ALPHA_SAT_FASTA):
            references.append(("alpha satellite", ALPHA_SAT_FASTA))

        if os.path.isfile(IMAGING_RISK_FASTA):
            references.append(("imaging-risk RNAs (tRNA/snRNA/7SL/7SK/mt-rRNA)", IMAGING_RISK_FASTA))

        if not references:
            self.logger.warning("No rDNA/satellite reference files found — skipping filter")
            Bio.SeqIO.write(sequences, self.output().path, format="fasta")
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            for ref_name, ref_fasta in references:
                db_path = os.path.join(tmpdir, os.path.basename(ref_fasta))
                blast_out = os.path.join(tmpdir, f"{os.path.basename(ref_fasta)}_blast.tsv")

                if not self._build_blast_db(ref_fasta, db_path):
                    continue

                with spinner(f"Screening probes against {ref_name}..."):
                    self._blast_probes(probe_fasta, db_path, blast_out, threads)

                dangerous = self._find_dangerous_probes(blast_out)
                if dangerous:
                    self.logger.debug(
                        f"Found {len(dangerous)} probes with dangerous near-matches "
                        f"to {ref_name}"
                    )
                all_dangerous |= dangerous

        if all_dangerous:
            self.logger.debug(
                f"Removing {len(all_dangerous)} probes with rDNA/satellite off-targets: "
                f"{sorted(all_dangerous)[:5]}{'...' if len(all_dangerous) > 5 else ''}"
            )

        candidates = [seq for seq in sequences if seq.id not in all_dangerous]
        util.log_and_check_candidates(
            self.logger, "FilterRibosomalRNA", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output().path, format="fasta")
