"""Off-target filters for alignments based on the reference genome sequence."""

import bisect
import logging
import os
import shutil
import subprocess

import Bio.Seq
import pandas as pd
import pysam

from . import constants, util
from .basic_filtering import compute_off_target_tm
from .config import GeneralConfig, ProbeConfig


class GenomeOffTargetFilterMixin:
    """Filter alignment hits by repeat masking and off-target duplex stability.

    Mixed into AlignProbeCandidates, which provides the filter settings.
    """

    logger = logging.getLogger("custom-logger")
    is_endogenous: bool
    max_off_targets: int
    off_target_min_tm: float

    def _run_dustmasker(self) -> str:
        """Run dustmasker on the reference genome and cache results."""
        genome = os.path.abspath(GeneralConfig().reference_genome)
        output_path = os.path.join(
            util.get_output_dir(),
            f"{util.get_gene_name()}_dustmasker.intervals",
        )
        if os.path.isfile(output_path):
            self.logger.debug(f"Using cached dustmasker output: {output_path}")
            return output_path

        from .console import spinner

        with spinner("Identifying repetitive regions (dustmasker)..."):
            subprocess.check_call(
                [
                    "dustmasker",
                    "-in", genome,
                    "-outfmt", "interval",
                    "-out", output_path,
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        return output_path

    @staticmethod
    def _load_masked_intervals(path: str) -> dict:
        """Parse dustmasker interval output into sorted interval lookup.

        Returns dict of {chrom: (sorted_starts, sorted_ends)}.
        """
        masked = {}
        current_chrom = None
        starts, ends = [], []

        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_chrom and starts:
                        masked[current_chrom] = (starts, ends)
                    current_chrom = line[1:].strip()
                    starts, ends = [], []
                elif " - " in line:
                    parts = line.split(" - ")
                    starts.append(int(parts[0]))
                    ends.append(int(parts[1]))
        if current_chrom and starts:
            masked[current_chrom] = (starts, ends)

        return masked

    def _is_in_masked_region(self, chrom: str, pos: int, length: int) -> bool:
        """Check if an alignment position overlaps a masked region."""
        if chrom not in self._masked_intervals:
            return False
        starts, ends = self._masked_intervals[chrom]
        # Find first interval that could overlap [pos, pos+length)
        idx = bisect.bisect_right(starts, pos) - 1
        for i in range(max(0, idx), min(len(starts), idx + 2)):
            if starts[i] < pos + length and ends[i] >= pos:
                return True
        return False

    def filter_masked_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove off-target hits in repetitive/masked regions before counting.

        This modifies the SAM DataFrame to exclude hits in masked regions,
        so the subsequent count-based filter sees fewer off-targets per probe.
        """
        if not shutil.which("dustmasker"):
            self.logger.warning(
                "dustmasker not found — skipping repeat masking. "
                "Install BLAST+ to enable."
            )
            return df_sam

        mask_path = self._run_dustmasker()
        self._masked_intervals = self._load_masked_intervals(mask_path)

        total_before = len(df_sam)
        is_masked = df_sam.apply(
            lambda row: self._is_in_masked_region(
                row["rname"], int(row["pos"]), len(str(row["seq"]))
            ),
            axis=1,
        )
        df_sam = df_sam[~is_masked]
        removed = total_before - len(df_sam)
        self.logger.debug(
            f"Removed {removed} alignments in masked/repetitive regions"
        )
        return df_sam

    def _is_significant_off_target(
        self,
        genome: pysam.FastaFile,
        hit: pd.Series,
        probe_seq: str,
        na: float,
        formamide: float,
    ) -> bool:
        """Check if an off-target hit is stable enough to count (Tm >= off_target_min_tm).

        Hits whose genomic target sequence cannot be fetched count as significant.
        """
        chrom = hit["rname"]
        pos = int(hit["pos"]) - 1  # SAM is 1-based, pysam is 0-based
        hit_len = len(hit["seq"])
        is_reverse = hit["flag"] == constants.SAM_FLAG_REVERSE

        try:
            target_seq = genome.fetch(chrom, pos, pos + hit_len).upper()
        except (ValueError, KeyError):
            return True

        if is_reverse:
            target_seq = str(Bio.Seq.Seq(target_seq).reverse_complement())

        tm = compute_off_target_tm(probe_seq, target_seq, na, formamide)
        return tm >= self.off_target_min_tm

    def filter_thermodynamic_off_targets(
        self, df_sam: pd.DataFrame, df_pre_count: pd.DataFrame
    ) -> pd.DataFrame:
        """Rescue probes whose off-targets are thermodynamically unstable.

        Re-examines probes that were removed by count-based filtering.
        For each off-target alignment, extracts the genomic target sequence
        and computes the predicted Tm. If all off-target Tms are below
        off_target_min_tm, the probe is rescued.

        Uses df_pre_count (the DataFrame after mask/intergenic filtering but
        before count-based filtering) so rescue evaluation respects the same
        pre-filters as the count-based step.
        """
        genome_path = os.path.abspath(GeneralConfig().reference_genome)

        # Index the genome if needed
        fai_path = genome_path + ".fai"
        if not os.path.isfile(fai_path):
            pysam.faidx(genome_path)

        genome = pysam.FastaFile(genome_path)

        # Find probes that were removed by count-based filtering
        kept_probes = set(df_sam["qname"].unique())
        all_probes = set(df_pre_count["qname"].unique())
        removed_probes = all_probes - kept_probes

        if not removed_probes:
            genome.close()
            return df_sam

        self.logger.debug(
            f"Thermodynamic re-evaluation of {len(removed_probes)} removed probes"
        )

        na = ProbeConfig().na_concentration
        formamide = ProbeConfig().formamide_concentration
        max_allowed = self.max_off_targets + (1 if self.is_endogenous else 0)
        rescued = []

        for probe_name in removed_probes:
            hits = df_pre_count[df_pre_count["qname"] == probe_name]
            probe_seq = hits.iloc[0]["seq"]

            # Evaluate each off-target hit (skip self-hit for endogenous)
            significant_count = sum(
                self._is_significant_off_target(genome, hit, probe_seq, na, formamide)
                for _, hit in hits.iterrows()
            )
            if significant_count <= max_allowed:
                rescued.append(probe_name)

        genome.close()

        if rescued:
            self.logger.debug(
                f"Rescued {len(rescued)} probes with thermodynamically "
                f"unstable off-targets (Tm < {self.off_target_min_tm}°C)"
            )
            # Add rescued probes back from the pre-count data (respects prior filters)
            rescued_rows = df_pre_count[df_pre_count["qname"].isin(rescued)]
            df_sam = pd.concat([df_sam, rescued_rows], ignore_index=True)

        return df_sam
