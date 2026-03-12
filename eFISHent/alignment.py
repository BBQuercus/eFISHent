"""Filter probes based on alignment score and uniqueness."""

import bisect
import logging
import os
import shutil
import subprocess

import Bio.Seq
import Bio.SeqIO
import luigi
import pandas as pd
import pysam

from . import constants
from . import util
from .basic_filtering import BasicFiltering
from .basic_filtering import compute_off_target_tm
from .config import GeneralConfig
from .config import ProbeConfig
from .config import SequenceConfig
from .indexing import BuildBowtieIndex
from .indexing import BuildBowtie2Index


class AlignProbeCandidates(luigi.Task):
    """Align all probe candidates to the reference genome and filter."""

    logger = logging.getLogger("custom-logger")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_endogenous = SequenceConfig().is_endogenous
        self.max_off_targets = ProbeConfig().max_off_targets
        self.no_alternative_loci = ProbeConfig().no_alternative_loci
        self.aligner = ProbeConfig().aligner
        self.encode_count_table = ProbeConfig().encode_count_table
        self.off_target_min_tm = ProbeConfig().off_target_min_tm
        self.mask_repeats = ProbeConfig().mask_repeats
        self.intergenic_off_targets = ProbeConfig().intergenic_off_targets
        self.filter_rrna = ProbeConfig().filter_rrna
        self.has_expression_filter = bool(
            self.encode_count_table
            and GeneralConfig().reference_annotation
            and self.is_endogenous
        )
        self.has_annotation = bool(GeneralConfig().reference_annotation)
        self.has_intergenic_filter = bool(
            self.intergenic_off_targets
            and self.has_annotation
            and self.is_endogenous
        )
        self.has_rrna_filter = bool(
            self.filter_rrna
            and self.has_annotation
            and self.is_endogenous
        )

    def requires(self):
        if self.aligner == "bowtie2":
            index_task = BuildBowtie2Index()
        else:
            index_task = BuildBowtieIndex()
        tasks = {"probes": BasicFiltering(), "bowtie": index_task}
        if self.has_expression_filter or self.has_intergenic_filter or self.has_rrna_filter:
            from .indexing import PrepareAnnotationFile
            tasks["annotation"] = PrepareAnnotationFile()
        return tasks

    def output(self):
        tasks = {
            "fasta": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_aligned.fasta"
                )
            ),
            "table": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(), f"{util.get_gene_name()}_counts.csv"
                )
            ),
        }
        return tasks

    @property
    def count_table(self) -> pd.DataFrame:
        """Read and verify a RNAseq normalized count table."""
        extension = os.path.splitext(self.fname_count)[1].lower()
        if extension not in (".tsv", ".csv", ".txt"):
            raise ValueError(
                "The count table must be provided as TSV, CSV, or TXT file. "
                f"Only found - {self.fname_count}"
            )
        sep = "\t" if extension in (".tsv", ".txt") else ","
        df = pd.read_csv(self.fname_count, sep=sep, index_col=0).reset_index().dropna()
        if not len(df.columns) >= 2:
            raise ValueError(
                "The count table must contain at least 2 columns."
                f"Only found - {df.columns}"
            )
        self.logger.debug(f"Mapping original the first two columns of - {df.columns}.")
        self.logger.debug(f"Read count table with {len(df)} entries.")
        self.logger.debug(f"Found dtypes - {df.dtypes}")
        self.logger.debug(f"First column - {df.iloc[0]}")

        df = df.iloc[:, :2]
        df.columns = ["gene_id", "count"]
        try:
            df = df.astype({"gene_id": str, "count": float})
        except ValueError:
            raise ValueError(
                "Could not convert count table columns to the right type. "
                "Please ensure the first two colums are of type string and number/float/int, respectively. "
            )

        # Remove potential version numbers
        df["gene_id"] = df["gene_id"].apply(lambda x: x.split(".")[0])
        return df

    @property
    def gtf_table(self) -> pd.DataFrame:
        """Parse prepared parquet-formatted gtf file."""
        return pd.read_parquet(self.fname_gtf)[constants.GTF_COLUMNS]

    @property
    def norm_table(self) -> pd.DataFrame:
        """Merge gtf and count table."""
        df = pd.merge(
            self.count_table, self.gtf_table, how="left", on="gene_id",
            validate="many_to_many",
        ).dropna()
        return df

    def align_probes(self, threads: int):
        """Align probes to the reference genome.

        Takes a fasta file as input (which will be converted to fastq).
        Saves a sam file with alignments.
        """
        # Convert fasta to fastq - bowtie only allows uses fastq inputs
        fname_fastq = self.fname_fasta.rstrip("a") + "q"
        with open(self.fname_fasta, "r") as fasta, open(fname_fastq, "w") as fastq:
            for record in Bio.SeqIO.parse(fasta, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                Bio.SeqIO.write(record, fastq, format="fastq")
        self.logger.debug(f"Converted fasta to fastq - {fname_fastq}")

        # Actual alignment with -m to only return alignments with >m hits
        alignments = str(max(self.max_off_targets + 1, 1))
        endo_exo_param = [] if self.is_endogenous else ["-m", alignments]
        args_bowtie = [
            "bowtie",
            fname_fastq,
            "-x",
            self.fname_genome,
            "--threads",
            str(threads),
            *endo_exo_param,
            "--sam",
            self.fname_sam,
            "--all",
        ]
        self.logger.debug(f"Running bowtie with - {' '.join(args_bowtie)}")
        from .console import spinner

        with spinner("Aligning probes to genome..."):
            subprocess.check_call(
                args_bowtie, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )

    def align_probes_bowtie2(self, threads: int):
        """Align probes using Bowtie2 with OligoMiner/Tigerfish parameters.

        Parameters are from OligoMiner (Beliveau et al., 2018, MIT license)
        and validated by Tigerfish (Kuo et al., 2024). Tuned for short
        oligonucleotide probe alignment with mismatch tolerance.
        """
        # Bowtie2 misparses FASTA descriptions containing angle brackets
        # (e.g., BioPython's "<unknown description>"). Write a clean FASTA
        # with IDs only to avoid read name corruption.
        fname_clean = self.fname_fasta + ".bt2.fa"
        with open(fname_clean, "w") as out:
            for record in Bio.SeqIO.parse(self.fname_fasta, "fasta"):
                out.write(f">{record.id}\n{record.seq}\n")

        k_value = str(max(self.max_off_targets + 2, 100))
        args = [
            "bowtie2",
            "-f",  # FASTA input
            "-x", self.fname_genome,
            "-U", fname_clean,
            "--threads", str(threads),
            "--local",  # Local alignment (soft-clip ends)
            "-D", "20",  # 20 seed extension attempts
            "-R", "3",  # 3 re-seeding rounds
            "-N", "1",  # 1 mismatch allowed in seed
            "-L", "20",  # Seed length 20bp
            "-i", "C,4",  # Seed interval: constant, every 4 bases
            "--score-min", "G,1,4",  # Min score: log-based threshold
            "-k", k_value,  # Report up to k alignments per probe
            "-S", self.fname_sam,
        ]
        # Only suppress unaligned for endogenous (exogenous needs unmapped reads)
        if self.is_endogenous:
            args.insert(-2, "--no-unal")
        self.logger.debug(f"Running bowtie2 with - {' '.join(args)}")
        from .console import spinner

        with spinner("Aligning probes to genome (Bowtie2)..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )
        os.remove(fname_clean)

    @staticmethod
    def parse_raw_pysam(sam: str) -> pd.DataFrame:
        """Convert string based pysam output to a DataFrame."""
        # Parse tab and newline delimited pysam output
        data = [row.split("\t")[:10] for row in sam.split("\n")]
        # end = 14 if is_endogenous else 12
        columns = constants.SAMFILE_COLUMNS[:10]

        # If empty [['']]
        if data[0][0]:
            return pd.DataFrame(data, columns=columns).dropna()
        return pd.DataFrame(columns=columns)

    def filter_unique_probes(self) -> tuple:
        """Filter sam file with probes based on alignment score and uniqueness.

        Filtering explanations:
            * mapq 60: uniquely mapped read, regardless of number of mismatches / indels
            * flag 4: unmapped read
        For Bowtie2:
            * Endogenous: keep all mapped reads, filter by hit count
            * Exogenous: keep only unmapped (flag 4)

        Returns (df_filtered, df_pre_count) where df_pre_count is the DataFrame
        after mask/intergenic filtering but before count-based filtering (used by
        thermodynamic rescue to evaluate probes against the same filtered hits).
        """
        if self.aligner == "bowtie2":
            if self.is_endogenous:
                # Keep all mapped reads, filter by hit count below
                filtered_sam = pysam.view(self.fname_sam, "--exclude-flags", "4")  # type: ignore
            else:
                filtered_sam = pysam.view(self.fname_sam, "--require-flags", "4")  # type: ignore
        else:
            flags = ["--min-MQ", "60"] if self.is_endogenous else ["--require-flags", "4"]
            self.logger.debug(f"Running samtools with - {' '.join(flags)}")
            filtered_sam = pysam.view(self.fname_sam, *flags)  # type: ignore

        df = self.parse_raw_pysam(filtered_sam)

        if self.no_alternative_loci:
            df = df[~df["rname"].str.contains("_alt")]

        # Remove hits in repetitive/masked regions before counting
        if self.mask_repeats and self.is_endogenous:
            df = self.filter_masked_off_targets(df)

        # Remove hits in intergenic regions before counting
        if self.has_intergenic_filter:
            df = self.filter_intergenic_off_targets(df)

        # Snapshot before count-based filtering (for thermodynamic rescue)
        df_pre_count = df.copy()

        # Filter to remove off-target rates
        if self.is_endogenous:
            group_sizes = df.groupby("qname").size()
            df = df[
                df["qname"].isin(
                    group_sizes[group_sizes <= self.max_off_targets + 1].index
                )
            ]
        return df, df_pre_count

    def get_most_expressed_genes(self, percentage: float) -> pd.Series:
        """Get gene IDs expressed above percentile threshold."""
        threshold = self.norm_table["count"].quantile(1 - percentage / 100)
        return self.norm_table[self.norm_table["count"] >= threshold]["gene_id"]

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

    def filter_rrna_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove probes with off-target hits on ribosomal RNA genes.

        rRNA is ~80% of cellular RNA — even weak binding causes intense
        background. Checks gene_biotype/gene_type for rRNA and Mt_rRNA.
        Any probe with at least one rRNA hit is vetoed entirely.
        """
        self.fname_gtf = self.input()["annotation"].path
        gtf_full = pd.read_parquet(self.fname_gtf)

        # Find the biotype column (differs between Ensembl/GENCODE)
        biotype_col = None
        for col in ["gene_biotype", "gene_type"]:
            if col in gtf_full.columns:
                biotype_col = col
                break

        if biotype_col is None:
            self.logger.debug(
                "No gene_biotype/gene_type column in GTF — skipping rRNA filter"
            )
            return df_sam

        # Get rRNA gene regions
        rrna_types = {"rRNA", "Mt_rRNA", "rRNA_pseudogene"}
        rrna_genes = gtf_full[
            (gtf_full[biotype_col].isin(rrna_types))
            & (gtf_full["feature"] == "gene")
        ][["seqname", "start", "end"]].drop_duplicates()

        if rrna_genes.empty:
            self.logger.debug("No rRNA genes found in GTF annotation")
            return df_sam

        self.logger.debug(f"Found {len(rrna_genes)} rRNA gene regions in GTF")

        # Build sorted interval lookup per chromosome for efficient overlap check
        rrna_by_chrom = {}
        for _, gene in rrna_genes.iterrows():
            chrom = gene["seqname"]
            if chrom not in rrna_by_chrom:
                rrna_by_chrom[chrom] = ([], [])
            rrna_by_chrom[chrom][0].append(gene["start"])
            rrna_by_chrom[chrom][1].append(gene["end"])

        def _hits_rrna(row):
            chrom = row["rname"]
            if chrom not in rrna_by_chrom:
                return False
            pos = int(row["pos"])
            starts, ends = rrna_by_chrom[chrom]
            return any(s <= pos <= e for s, e in zip(starts, ends))

        hits_rrna = df_sam.apply(_hits_rrna, axis=1)
        probes_to_remove = set(df_sam.loc[hits_rrna, "qname"])

        if probes_to_remove:
            self.logger.debug(
                f"Removing {len(probes_to_remove)} probes hitting rRNA off-targets"
            )
            df_sam = df_sam[~df_sam["qname"].isin(probes_to_remove)]

        return df_sam

    def filter_intergenic_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove off-target hits in intergenic regions before counting.

        Hits not overlapping any annotated gene are unlikely to produce RNA
        that the FISH probe could bind, so they are excluded from off-target
        counts. Uses the prepared GTF annotation.
        """
        self.fname_gtf = self.input()["annotation"].path
        gtf = self.gtf_table
        genes = gtf[gtf["feature"] == "gene"]

        total_before = len(df_sam)

        def _overlaps_gene(row):
            chrom = row["rname"]
            pos = int(row["pos"])
            return (
                (genes["seqname"] == chrom)
                & (genes["start"] <= pos)
                & (genes["end"] >= pos)
            ).any()

        is_genic = df_sam.apply(_overlaps_gene, axis=1)
        df_sam = df_sam[is_genic]
        removed = total_before - len(df_sam)
        self.logger.debug(
            f"Removed {removed} alignments in intergenic regions"
        )
        return df_sam

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
        rescued = []

        for probe_name in removed_probes:
            hits = df_pre_count[df_pre_count["qname"] == probe_name]
            probe_seq = hits.iloc[0]["seq"]

            # Evaluate each off-target hit (skip self-hit for endogenous)
            significant_count = 0
            for _, hit in hits.iterrows():
                chrom = hit["rname"]
                pos = int(hit["pos"]) - 1  # SAM is 1-based, pysam is 0-based
                hit_len = len(hit["seq"])
                is_reverse = hit["flag"] == constants.SAM_FLAG_REVERSE

                try:
                    target_seq = genome.fetch(chrom, pos, pos + hit_len).upper()
                except (ValueError, KeyError):
                    significant_count += 1
                    continue

                if is_reverse:
                    target_seq = str(Bio.Seq.Seq(target_seq).reverse_complement())

                tm = compute_off_target_tm(probe_seq, target_seq, na, formamide)
                if tm >= self.off_target_min_tm:
                    significant_count += 1

            max_allowed = self.max_off_targets + (1 if self.is_endogenous else 0)
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

    def filter_expression_off_targets(self, df_sam: pd.DataFrame) -> pd.DataFrame:
        """Remove probes aligning to highly expressed off-target genes."""
        self.fname_count = self.encode_count_table
        self.fname_gtf = self.input()["annotation"].path

        expressed_genes = self.get_most_expressed_genes(
            ProbeConfig().max_expression_percentage
        )
        self.logger.debug(
            f"Found {len(expressed_genes)} highly expressed genes "
            f"(top {ProbeConfig().max_expression_percentage}%)"
        )

        # Join alignments with gene annotations by position overlap
        gtf = self.gtf_table
        gtf = gtf[gtf["feature"] == "gene"]

        # For each alignment, find overlapping genes
        probes_to_remove = set()
        for _, row in df_sam.iterrows():
            chrom = row["rname"]
            pos = int(row["pos"])
            # Find genes overlapping this position
            overlapping = gtf[
                (gtf["seqname"] == chrom)
                & (gtf["start"] <= pos)
                & (gtf["end"] >= pos)
            ]
            for _, gene in overlapping.iterrows():
                gene_id = gene["gene_id"]
                if gene_id in expressed_genes.values:
                    probes_to_remove.add(row["qname"])

        if probes_to_remove:
            self.logger.debug(
                f"Removing {len(probes_to_remove)} probes hitting highly expressed off-targets"
            )
            df_sam = df_sam[~df_sam["qname"].isin(probes_to_remove)]

        return df_sam

    def run(self):
        util.log_stage_start(self.logger, "AlignProbeCandidates")
        # Naming
        self.fname_fasta = self.input()["probes"].path
        self.fname_sam = os.path.splitext(self.fname_fasta)[0] + ".sam"
        self.fname_genome = util.get_genome_name()

        # Alignment
        if self.aligner == "bowtie2":
            self.align_probes_bowtie2(threads=GeneralConfig().threads)
        else:
            self.align_probes(threads=GeneralConfig().threads)

        # Filtering
        df_sam, df_pre_count = self.filter_unique_probes()

        # Thermodynamic off-target filtering (rescue probes with unstable off-targets)
        if self.off_target_min_tm > 0 and self.is_endogenous:
            df_sam = self.filter_thermodynamic_off_targets(df_sam, df_pre_count)

        # rRNA off-target filtering
        if self.has_rrna_filter:
            df_sam = self.filter_rrna_off_targets(df_sam)

        # Expression-weighted filtering
        if self.has_expression_filter:
            df_sam = self.filter_expression_off_targets(df_sam)

        df_sam.to_csv(self.output()["table"].path, index=False)

        # Save output
        sequences = list(Bio.SeqIO.parse(self.fname_fasta, "fasta"))
        candidates = [seq for seq in sequences if seq.id in df_sam["qname"].unique()]
        self.logger.debug(df_sam["qname"].values)
        util.log_and_check_candidates(
            self.logger, "AlignProbeCandidates", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, format="fasta")
