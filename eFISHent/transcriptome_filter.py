"""Filter probes based on BLAST hits against a reference transcriptome.

BLAST parameters adapted from TrueProbes (Neuert Lab, 2025). Uses blastn
(not megablast) with word_size 7 to detect short partial matches on probes.
The approach of BLASTing probes against the transcriptome (rather than genome)
catches off-target binding to expressed transcripts that alignment-only methods miss.
"""

import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import Bio.Seq
import Bio.SeqIO
import luigi
import numpy as np
import pandas as pd

from . import util
from .config import GeneralConfig, ProbeConfig

# Sensitivity settings shared by the filter and the final-probe report.
# With fewer reported targets, short perfect matches (e.g. 16 nt, e~1) outrank
# the 18+ nt one-mismatch hits (e~14) that the off-target rule counts and push
# them out of the result, so qualifying off-targets go unseen. Values follow
# TrueProbes (Neuert Lab, 2025): word 7, evalue 1000, 1000 alignments.
BLAST_EVALUE = "1000"
BLAST_MAX_TARGETS = "1000"


def antisense_hits(df_blast: pd.DataFrame) -> pd.DataFrame:
    """Keep hits where the probe is complementary to the transcript.

    The transcriptome FASTA is in transcript (sense) orientation and probes are
    reverse complements of their target, so a probe can only hybridize to a
    transcript it aligns to plus/minus (sstart > send). A plus/plus hit means
    the probe has the transcript's own sequence (e.g. an antisense RNA
    overlapping the target) and binds only its absent complement.
    """
    return df_blast[df_blast["sstart"] > df_blast["send"]]


# Hybridization temperature the off-target burden is scored against. Probe Tms
# are formamide-corrected, so this is the actual incubation temperature.
HYBRIDIZATION_TEMP = 37.0
# Shortest hit that counts towards the burden (TrueProbes ignores < 15 nt).
MIN_BURDEN_MATCH_LENGTH = 15


def hit_melting_temps(
    df: pd.DataFrame,
    probe_seqs: dict[str, str],
    na_concentration: float,
    formamide_concentration: float,
) -> np.ndarray:
    """Predict the duplex Tm of each BLAST hit.

    Tm of the matched probe segment (qstart..qend) as a perfect RNA/DNA
    hybrid, minus 2.5 °C per mismatch or gap (as compute_off_target_tm).
    Segment Tms are cached, since most hits share a few segments.
    """
    from .basic_filtering import get_melting_temp

    cache: dict[tuple[str, int, int], float] = {}
    tms = np.empty(len(df))
    rows = zip(df["qseqid"], df["qstart"], df["qend"], df["mismatch"], df["gapopen"])
    for i, (probe, qstart, qend, mismatch, gaps) in enumerate(rows):
        key = (probe, int(qstart), int(qend))
        if key not in cache:
            lo, hi = sorted((key[1], key[2]))
            segment = Bio.Seq.Seq(probe_seqs[probe][lo - 1 : hi])
            cache[key] = get_melting_temp(segment, na_concentration, formamide_concentration)
        tms[i] = cache[key] - 2.5 * (mismatch + gaps)
    return tms


def binding_weight(tm: np.ndarray, temperature: float = HYBRIDIZATION_TEMP) -> np.ndarray:
    """Logistic estimate (0-1) that a duplex with this Tm is bound at temperature."""
    return 1.0 / (1.0 + np.exp(-(np.asarray(tm, dtype=float) - temperature) / 2.0))


_BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
]

# Cross-hybridization: a contiguous (ungapped) match at >=95% identity that
# covers >=85% of the probe length. For a 20nt probe: >=17nt. For a 24nt
# probe: >=20nt.
CROSS_HYB_MIN_FRACTION = 0.85
CROSS_HYB_MIN_IDENTITY = 95.0


@dataclass(frozen=True)
class BlastSearch:
    """BLAST settings shared by all queries of one transcriptome search."""

    task: str
    db: str
    perc_identity: float

    def args(self, query: str, out: str, threads: int) -> list[str]:
        """Build BLAST command-line arguments."""
        return [
            "blastn",
            "-task", self.task,
            "-query", query,
            "-db", self.db,
            "-evalue", BLAST_EVALUE,
            "-word_size", "7",
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "1",
            "-penalty", "-3",
            "-dust", "no",
            "-max_target_seqs", BLAST_MAX_TARGETS,
            "-perc_identity", str(self.perc_identity),
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore",
            "-out", out,
        ]


def _run_blast_chunk(args):
    """Run a single BLAST process (top-level function for pickling)."""
    subprocess.check_call(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def _concatenate_files(paths: list[str], out: str) -> None:
    """Concatenate the non-empty files in paths into out."""
    with open(out, "w") as fout:
        for path in paths:
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                with open(path) as fin:
                    fout.write(fin.read())


def _read_blast_hits(blast_out: str) -> pd.DataFrame:
    """Read antisense BLAST hits with their effective match length.

    Effective match length = alignment length - gaps (per TrueProbes).
    """
    if os.path.getsize(blast_out) > 0:
        df_blast = pd.read_csv(blast_out, sep="\t", header=None, names=_BLAST_COLUMNS)
    else:
        df_blast = pd.DataFrame(columns=_BLAST_COLUMNS)
    df_blast = antisense_hits(df_blast)
    df_blast["effective_len"] = df_blast["length"] - df_blast["gapopen"]
    return df_blast


def _name_pattern(name: str) -> str:
    """Regex matching name as a delimited token of a subject ID."""
    return rf"(?:^|[|_\-.\s]){re.escape(name)}(?:$|[|_\-.\s])"


def _clean_gene_name(name: str) -> str:
    """Strip common suffixes and species codes from a gene name.

    E.g., "eif2b1_cds" -> "eif2b1", "mettl3_cds_odd" -> "mettl3",
    "chd4_mm" -> "chd4", "taf10_hs" -> "taf10".
    """
    name = re.sub(r"[_-]?(cds|mrna|rna|transcript|seq|odd|even).*$", "", name)
    return re.sub(r"[_-](?:mm|hs|rn|dm|dr|ce|sc)$", "", name)


class SelfHitMatcher:
    """Decide whether a BLAST subject is a transcript of the target gene.

    A subject is a self-hit if one of the gene name patterns matches its ID,
    or if its ID (with or without version suffix) is a known target transcript.
    """

    def __init__(
        self, patterns: list[str], transcript_ids: set[str], cli_gene_name: str = ""
    ):
        self.patterns = patterns
        self.transcript_ids = transcript_ids
        self.cli_gene_name = cli_gene_name
        self.match_cli_substring = False

    def promote(self, transcript_ids: set[str]) -> None:
        """Count transcripts as self, and --gene-name anywhere in an ID."""
        self.transcript_ids.update(transcript_ids)
        self.match_cli_substring = bool(self.cli_gene_name)

    def _matches_name(self, sid_lower: str) -> bool:
        if any(re.search(p, sid_lower) for p in self.patterns):
            return True
        return self.match_cli_substring and self.cli_gene_name in sid_lower

    def __call__(self, sseqid: str) -> bool:
        sid_lower = sseqid.lower()
        if self._matches_name(sid_lower):
            return True
        base = sid_lower.split(".")[0]
        return sid_lower in self.transcript_ids or base in self.transcript_ids


def _phantom_self_transcripts(
    df_blast: pd.DataFrame, n_probes: int, is_self_hit: SelfHitMatcher
) -> set[str]:
    """Find non-self transcripts that are alternative annotations of the target.

    If a non-self transcript hits >=50% of all probes at 100% identity it is
    almost certainly a read-through or pseudogene overlay of the target, not a
    true off-target.
    """
    if n_probes == 0:
        return set()
    perfect_hits = df_blast[(df_blast["pident"] >= 100.0) & (df_blast["gapopen"] == 0)]
    if perfect_hits.empty:
        return set()
    phantom_counts = (
        perfect_hits[~perfect_hits["sseqid"].apply(is_self_hit)]
        .groupby("sseqid")["qseqid"]
        .nunique()
    )
    return set(phantom_counts[phantom_counts >= n_probes * 0.5].index.str.lower())


def _off_target_hits(
    df_blast: pd.DataFrame,
    is_self_hit: SelfHitMatcher,
    probe_seqs: dict[str, str],
    config: ProbeConfig,
) -> pd.DataFrame:
    """Off-target hits >= 15 nt with their predicted duplex Tm.

    Self-hit status is resolved once per subject, not per BLAST row.
    """
    self_subjects = {sid: is_self_hit(sid) for sid in df_blast["sseqid"].unique()}
    df_ot = df_blast[
        (df_blast["effective_len"] >= MIN_BURDEN_MATCH_LENGTH)
        & ~df_blast["sseqid"].map(self_subjects).astype(bool)
    ].copy()
    df_ot["tm"] = hit_melting_temps(
        df_ot, probe_seqs, config.na_concentration, config.formamide_concentration
    )
    return df_ot


def _count_off_targets(
    df_hits: pd.DataFrame, is_self_hit: SelfHitMatcher
) -> dict[str, int]:
    """Count the unique off-target transcripts hit by each probe."""
    off_target_counts = {}
    for probe_id, group in df_hits.groupby("qseqid"):
        off_targets = group[~group["sseqid"].apply(is_self_hit)]
        off_target_counts[probe_id] = len(off_targets["sseqid"].unique())
    return off_target_counts


def _cross_hyb_threshold(probe_len: int, override: int) -> int:
    """Shortest contiguous match that counts as cross-hybridization."""
    if override > 0:
        return override
    return max(15, int(CROSS_HYB_MIN_FRACTION * probe_len))


def _cross_hybridizing_probes(
    df_blast: pd.DataFrame,
    probe_lengths: dict[str, int],
    is_self_hit: SelfHitMatcher,
    config: ProbeConfig,
) -> set[str]:
    """Find probes with a long contiguous near-perfect match to an off-target.

    These near-perfect matches to abundant or spatially concentrated RNAs
    create noisy images.
    """
    df_contiguous = df_blast[
        (df_blast["pident"] >= CROSS_HYB_MIN_IDENTITY)
        & (df_blast["gapopen"] == 0)
    ]
    cross_hyb_probes = set()
    for probe_id, group in df_contiguous.groupby("qseqid"):
        threshold = _cross_hyb_threshold(
            probe_lengths.get(probe_id, config.min_length), config.cross_hyb_min_length
        )
        hits = group[group["effective_len"] >= threshold]
        if not hits.empty and not hits[~hits["sseqid"].apply(is_self_hit)].empty:
            cross_hyb_probes.add(probe_id)
    return cross_hyb_probes


class BuildTranscriptomeBlastDB(luigi.Task):
    """Build a BLAST nucleotide database from a transcriptome FASTA."""

    logger = logging.getLogger("custom-logger")

    def output(self):
        txome = os.path.abspath(GeneralConfig().reference_transcriptome)
        return luigi.LocalTarget(txome + ".nsq")

    def run(self):
        from .console import spinner

        txome = os.path.abspath(GeneralConfig().reference_transcriptome)
        args = [
            "makeblastdb",
            "-in", txome,
            "-dbtype", "nucl",
            "-out", txome,
        ]
        self.logger.debug(f"Running makeblastdb with - {' '.join(args)}")
        with spinner("Building BLAST database..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )
        self.logger.info("Finished building BLAST transcriptome database.")


class TranscriptomeFiltering(luigi.Task):
    """Filter probes with significant off-target hits in the transcriptome."""

    logger = logging.getLogger("custom-logger")

    def requires(self):
        from .secondary_structure import SecondaryStructureFiltering
        return {
            "probes": SecondaryStructureFiltering(),
            "blastdb": BuildTranscriptomeBlastDB(),
        }

    def output(self):
        return {
            "fasta": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(),
                    f"{util.get_gene_name()}_txome.fasta",
                )
            ),
            "table": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(),
                    f"{util.get_gene_name()}_txome_hits.csv",
                )
            ),
            "burden": luigi.LocalTarget(
                os.path.join(
                    util.get_output_dir(),
                    f"{util.get_gene_name()}_txome_burden.csv",
                )
            ),
        }

    def run(self):
        util.log_stage_start(self.logger, "TranscriptomeFiltering")
        probes_input = self.input()["probes"]
        # SecondaryStructureFiltering returns a single LocalTarget
        probe_fasta = probes_input.path if hasattr(probes_input, "path") else probes_input["fasta"].path
        sequences = list(Bio.SeqIO.parse(probe_fasta, "fasta"))
        config = ProbeConfig()

        blast_out = os.path.splitext(probe_fasta)[0] + "_blast.tsv"
        self._run_blast(probe_fasta, sequences, blast_out, config)
        df_blast = _read_blast_hits(blast_out)
        df_hits = self._counted_hits(df_blast, config)
        is_self_hit, tx_map = self._build_self_hit_matcher(df_blast, len(sequences))

        probe_seqs = {seq.id: str(seq.seq) for seq in sequences}
        df_ot = _off_target_hits(df_blast, is_self_hit, probe_seqs, config)
        self._write_burden(df_ot, sequences, tx_map)
        df_hits = self._apply_tm_gate(df_hits, df_ot, config.off_target_min_tm)

        off_target_counts = _count_off_targets(df_hits, is_self_hit)
        max_off = config.max_transcriptome_off_targets
        candidates = [
            seq for seq in sequences
            if off_target_counts.get(seq.id, 0) <= max_off
        ]

        probe_lengths = {seq.id: len(seq.seq) for seq in sequences}
        cross_hyb_probes = _cross_hybridizing_probes(
            df_blast, probe_lengths, is_self_hit, config
        )
        candidates = self._handle_cross_hybridization(
            candidates, cross_hyb_probes, config.reject_cross_hybridization
        )

        # Save hits table
        df_hits.to_csv(self.output()["table"].path, index=False)

        util.log_and_check_candidates(
            self.logger, "TranscriptomeFiltering", len(candidates), len(sequences)
        )
        Bio.SeqIO.write(candidates, self.output()["fasta"].path, format="fasta")

        # Clean up temp blast output
        if os.path.exists(blast_out):
            os.remove(blast_out)

    def _run_blast(self, probe_fasta, sequences, blast_out, config):
        """BLAST all probes against the transcriptome into blast_out."""
        txome_db = os.path.abspath(GeneralConfig().reference_transcriptome)
        self.logger.info(f"BLAST against {os.path.basename(txome_db)}")

        max_probe_len = config.max_length
        search = BlastSearch(
            # Use blastn-short for probe-length queries (<30bp) — optimized
            # seed strategy and scoring for short sequences
            task="blastn-short" if max_probe_len <= 30 else "blastn",
            db=txome_db,
            perc_identity=max(config.blast_identity_threshold, 100.0 * 15 / max_probe_len),
        )
        threads = GeneralConfig().threads

        # Parallel BLAST: split probes into chunks and run independent BLAST
        # processes. This scales better than BLAST's own -num_threads because
        # each process gets its own I/O and memory pipeline.
        n_chunks = max(1, min(threads, math.ceil(len(sequences) / 500)))
        if n_chunks <= 1:
            self._blast_single(search, probe_fasta, blast_out, threads)
        else:
            self._blast_chunked(search, sequences, blast_out, n_chunks)

    def _blast_single(self, search, probe_fasta, blast_out, threads):
        """BLAST a small probe set in a single call."""
        from .console import spinner

        args = search.args(probe_fasta, blast_out, threads)
        self.logger.debug(f"Running BLAST with - {' '.join(args)}")
        with spinner("BLASTing probes against transcriptome..."):
            subprocess.check_call(
                args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
            )

    def _blast_chunked(self, search, sequences, blast_out, n_chunks):
        """Split probes into chunks and BLAST them in parallel processes."""
        from .console import spinner

        self.logger.debug(
            f"Splitting {len(sequences)} probes into {n_chunks} chunks "
            f"for parallel BLAST"
        )
        chunk_size = math.ceil(len(sequences) / n_chunks)
        tmpdir = tempfile.mkdtemp(prefix="efishent_blast_")
        chunk_outs = []
        all_args = []

        try:
            for i in range(n_chunks):
                chunk = sequences[i * chunk_size : (i + 1) * chunk_size]
                if not chunk:
                    continue
                chunk_fa = os.path.join(tmpdir, f"chunk_{i}.fasta")
                chunk_out = os.path.join(tmpdir, f"chunk_{i}.tsv")
                Bio.SeqIO.write(chunk, chunk_fa, "fasta")
                chunk_outs.append(chunk_out)
                # Each chunk gets 1 thread — parallelism comes from processes
                all_args.append(search.args(chunk_fa, chunk_out, 1))

            with spinner(
                f"BLASTing probes against transcriptome ({n_chunks} parallel)..."
            ):
                with ProcessPoolExecutor(max_workers=n_chunks) as pool:
                    list(pool.map(_run_blast_chunk, all_args))

            _concatenate_files(chunk_outs, blast_out)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def _counted_hits(self, df_blast, config):
        """Hits long and similar enough to count as off-targets.

        Default min_match_len: max(18, 0.8 * min_probe_length), per TrueProbes.
        """
        min_match_len = config.min_blast_match_length
        if min_match_len == 0:
            min_match_len = max(18, int(0.8 * config.min_length))
        self.logger.debug(f"Using min_blast_match_length = {min_match_len}")
        return df_blast[
            (df_blast["effective_len"] >= min_match_len)
            & (df_blast["pident"] >= config.blast_identity_threshold)
        ]

    def _build_self_hit_matcher(self, df_blast, n_probes):
        """Build the self-hit predicate and the transcript-to-gene mapping.

        Strategy 1: gene name in subject ID (works for gffread-style headers).
        Strategy 2: GTF-based transcript-to-gene mapping (works for Ensembl IDs).
        """
        from .config import SequenceConfig

        gene_name_raw = util.get_gene_name(hashed=False).lower()
        gene_name_clean = _clean_gene_name(gene_name_raw)
        # Include the raw --gene-name CLI value (e.g., "URA3") which may differ
        # from the file-based gene_name_raw (e.g., "saccharomyces_cerevisiae_ura3")
        cli_gene_name = (SequenceConfig().gene_name or "").strip().lower()

        search_names = {gene_name_raw, gene_name_clean}
        if cli_gene_name:
            search_names.add(cli_gene_name)
        tx_map, self_transcript_ids = self._gtf_self_transcripts(
            search_names, gene_name_clean
        )

        # The CLI gene name matches e.g. "URA3" in "YEL021W_URA3_mRNA"
        names = [gene_name_raw]
        if gene_name_clean != gene_name_raw:
            names.append(gene_name_clean)
        if cli_gene_name:
            names.append(cli_gene_name)
        is_self_hit = SelfHitMatcher(
            [_name_pattern(name) for name in names], self_transcript_ids, cli_gene_name
        )

        phantom_self = _phantom_self_transcripts(df_blast, n_probes, is_self_hit)
        if phantom_self:
            self.logger.debug(
                f"Auto-promoted {len(phantom_self)} transcripts to self-hit "
                f"(≥50% of probes hit at 100% identity): {phantom_self}"
            )
            is_self_hit.promote(phantom_self)
        return is_self_hit, tx_map

    def _gtf_self_transcripts(self, search_names, gene_name_clean):
        """Map transcripts to genes with the GTF and find the target's transcripts.

        Returns the transcript-to-gene mapping and the lowercase IDs of all
        transcripts whose gene is one of search_names (both empty without GTF).
        """
        self_transcript_ids: set[str] = set()
        tx_map: dict[str, str] = {}
        gtf_path = GeneralConfig().reference_annotation
        if not gtf_path:
            return tx_map, self_transcript_ids
        try:
            from .gene_annotation import build_transcript_gene_map

            tx_map = build_transcript_gene_map(gtf_path)
            self_transcript_ids.update(
                tid.lower() for tid, gname in tx_map.items()
                if gname.lower() in search_names
            )
            if self_transcript_ids:
                self.logger.debug(
                    f"GTF mapping: found {len(self_transcript_ids)} transcripts "
                    f"for target gene '{gene_name_clean}'"
                )
        except Exception as e:
            self.logger.debug(f"GTF mapping failed, using name-based exclusion: {e}")
        return tx_map, self_transcript_ids

    def _write_burden(self, df_ot, sequences, tx_map):
        """Write each probe's off-target burden.

        Per off-target gene, the likeliest-bound hit counts; summed over genes.
        """
        df_ot["gene"] = [
            tx_map.get(sid, tx_map.get(sid.split(".")[0], sid)) for sid in df_ot["sseqid"]
        ]
        df_ot["weight"] = binding_weight(df_ot["tm"].values)
        burden = (
            df_ot.groupby(["qseqid", "gene"])["weight"].max().groupby("qseqid").sum()
        )
        pd.DataFrame(
            {"name": [seq.id for seq in sequences]}
        ).assign(
            burden=lambda d: d["name"].map(burden).fillna(0.0).round(4)
        ).to_csv(self.output()["burden"].path, index=False)

    def _apply_tm_gate(self, df_hits, df_ot, min_tm):
        """Drop counted hits whose duplex Tm is below min_tm.

        With --off-target-min-tm, hits that cannot form a stable duplex don't
        count towards max_transcriptome_off_targets.
        """
        if min_tm <= 0:
            return df_hits
        weak = df_ot.index[df_ot["tm"] < min_tm]
        if len(weak):
            n_before = df_hits["qseqid"].nunique()
            df_hits = df_hits.drop(index=df_hits.index.intersection(weak))
            self.logger.debug(
                f"Ignored off-target hits with Tm < {min_tm} °C "
                f"({n_before} -> {df_hits['qseqid'].nunique()} probes with counted hits)"
            )
        return df_hits

    def _handle_cross_hybridization(self, candidates, cross_hyb_probes, reject):
        """Reject (or warn about) probes with strong cross-hybridization."""
        if not cross_hyb_probes:
            return candidates
        if reject:
            self.logger.debug(
                f"Rejecting {len(cross_hyb_probes)} probes with strong "
                f"cross-hybridization (>=85% of probe length contiguous at >=95% identity)"
            )
            return [seq for seq in candidates if seq.id not in cross_hyb_probes]
        self.logger.debug(
            f"Cross-hybridization warning: {len(cross_hyb_probes)} probes have "
            f">=85% contiguous match at >=95% identity to off-target transcripts"
        )
        return candidates
