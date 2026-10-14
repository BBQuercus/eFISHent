"""Non-fatal warnings about risky parameter combinations."""

import argparse
import logging
import os
import sys


def _probe_parameter_warnings(args: argparse.Namespace) -> list[str]:
    """Warn about probe design thresholds that may yield few probes."""
    warnings = []

    # Narrow TM window
    tm_window = args.max_tm - args.min_tm
    if tm_window < 10:
        warnings.append(
            f"TM window is very narrow ({tm_window:.0f}\u00b0C) \u2014 may yield few probes.\n"
            "  Consider widening to at least 10\u00b0C."
        )

    # High formamide with low TM range
    if args.formamide_concentration > 30 and args.max_tm < 60:
        warnings.append(
            f"Formamide {args.formamide_concentration:.0f}% with max TM {args.max_tm:.0f}\u00b0C "
            "may filter all probes.\n"
            "  Consider --min-tm 50 --max-tm 70 for high formamide."
        )

    # K-mer length close to probe length
    if args.kmer_length >= args.min_length - 3:
        warnings.append(
            f"K-mer length ({args.kmer_length}) is close to min probe length ({args.min_length}).\n"
            f"  This may cause excessive filtering. Consider --kmer-length {max(10, args.min_length - 6)}."
        )

    # Single very short probe length
    if args.min_length == args.max_length and args.max_length < 20:
        warnings.append(
            f"Single probe length of {args.min_length}nt is very short.\n"
            "  Consider allowing a range (e.g., --min-length 20 --max-length 25)."
        )

    # Narrow GC window
    gc_window = args.max_gc - args.min_gc
    if gc_window < 20:
        warnings.append(
            f"GC window is very narrow ({gc_window:.0f}%) \u2014 may yield few probes.\n"
            "  Consider widening to at least 20%."
        )

    # Very strict secondary structure filter
    if args.max_deltag > -3:
        warnings.append(
            f"deltaG threshold ({args.max_deltag} kcal/mol) is very strict.\n"
            "  Consider --max-deltag -10 for a more permissive filter."
        )

    # Large spacing
    if args.spacing > 10:
        warnings.append(
            f"Probe spacing of {args.spacing}nt is large \u2014 may reduce coverage.\n"
            "  Consider --spacing 2 for denser probe tiling."
        )

    return warnings


def _prebuilt_genome_warnings(args: argparse.Namespace) -> list[str]:
    """Warn about slow paths when using pre-built genome indices."""
    if not getattr(args, "genome", None):
        return []

    warnings = []

    # Pre-built genome + non-default kmer_length: forced jellyfish rebuild
    # genome_15.jf is the only pre-built jellyfish index; any other kl builds
    # from the full genome and takes 10-30+ min on large genomes (longer on NFS).
    if args.kmer_length != 15:
        warnings.append(
            f"Using --genome with -kl {args.kmer_length} (not 15): "
            "the pre-built package only includes a 15-mer jellyfish index.\n"
            f"  eFISHent will check HuggingFace for a pre-built {args.kmer_length}-mer index "
            "and download it automatically if available.\n"
            "  If not available, the index will be built locally from the full genome "
            "(slow on large genomes / NFS — consider keeping -kl 15)."
        )

    # Pre-built genome on NFS: bowtie2 reads ~6 GB index; suggest cache on scratch.
    # Only shown on Linux (where HPC/NFS setups are the norm; macOS home is local).
    if not getattr(args, "index_cache_dir", None) and sys.platform.startswith("linux"):
        warnings.append(
            "Using --genome: indices default to ~/.local/efishent/indices/.\n"
            "  On HPC systems where HOME is NFS-mounted, genome alignment can take "
            "30-60+ min instead of ~1 min.\n"
            "  Add --index-cache-dir /path/to/scratch to store indices on a fast "
            "local filesystem."
        )

    return warnings


def _long_probe_kmer_warning(args: argparse.Namespace) -> str | None:
    """Warn when long probes meet a strict k-mer filter.

    A 50-mer has ~36 sub-15-mers, so the probability of any one exceeding
    max_kmers=5 is much higher than for a 20-24-mer (the smfish preset
    default target length).
    """
    if args.min_length <= 30 or args.max_kmers > 5:
        return None
    suggested = max(10, round(args.max_kmers * args.min_length / 24))
    return (
        f"Long probes ({args.min_length}–{args.max_length} nt) with strict k-mer filter "
        f"(--max-kmers {args.max_kmers}) may remove most probes.\n"
        f"  A {args.min_length}-mer contains ~{args.min_length - 14} sub-15-mers; "
        "the chance of any one exceeding the count threshold is much higher than for short probes.\n"
        f"  Consider --max-kmers {suggested} for long probe designs."
    )


def _sequence_gc_percent(seq_file: str) -> float:
    """Return the GC content (%) of the first record in a FASTA file."""
    import Bio.SeqIO
    import Bio.SeqUtils

    seq = next(Bio.SeqIO.parse(seq_file, "fasta"))
    return float(Bio.SeqUtils.gc_fraction(seq.seq)) * 100


def _high_gc_target_warning(args: argparse.Namespace) -> str | None:
    """Warn about high-GC targets, computing gene GC from the sequence file if available."""
    seq_file = getattr(args, "sequence_file", None)
    if not seq_file or not os.path.isfile(seq_file):
        return None

    try:
        gene_gc = _sequence_gc_percent(seq_file)
    except Exception as e:
        logging.getLogger("custom-logger").debug(
            f"Skipping high-GC target check for {seq_file}: {e}"
        )
        return None

    if gene_gc <= 55:
        return None
    return (
        f"High-GC target detected ({gene_gc:.1f}% GC) \u2014 elevated risk of "
        "nuclear/nucleolar background.\n"
        "  Probes from high-GC genes cross-hybridize with rRNA (~60-67% GC).\n"
        "  Consider tightening --max-gc (e.g., 60) or using --max-cpg-fraction 0.10."
    )


def validate_parameter_warnings(args: argparse.Namespace) -> list[str]:
    """Check for risky parameter combinations and return warning messages.

    These are non-fatal warnings — the pipeline will still run.
    """
    # Skip for non-pipeline modes
    if args.build_indices or args.analyze_probeset:
        return []

    warnings = _probe_parameter_warnings(args) + _prebuilt_genome_warnings(args)
    for warning in (_long_probe_kmer_warning(args), _high_gc_target_warning(args)):
        if warning:
            warnings.append(warning)
    return warnings
