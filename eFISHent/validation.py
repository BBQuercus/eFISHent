"""Validation of command line argument values and combinations."""

import argparse
import os


def string_to_bool(value) -> bool:
    """Workaround for using typed boolean values as arguments."""
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif value.lower() in ("no", "false", "f", "n", "0"):
        return False

    raise argparse.ArgumentTypeError(
        f"Invalid boolean value '{value}'. Use: true/false, yes/no, 1/0"
    )


def positive_int(value) -> int:
    """Validate that value is a positive integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'") from None
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"Value must be >= 1, got {ivalue}")
    return ivalue


def non_negative_int(value) -> int:
    """Validate that value is a non-negative integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'") from None
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"Value must be >= 0, got {ivalue}")
    return ivalue


def percentage(value) -> float:
    """Validate that value is a percentage (0-100)."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid number: '{value}'") from None
    if not 0 <= fvalue <= 100:
        raise argparse.ArgumentTypeError(f"Percentage must be 0-100, got {fvalue}")
    return fvalue


def existing_file(value) -> str:
    """Validate that file exists."""
    if not value:
        return value
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"File not found: '{value}'")
    return value


def existing_fasta_file(value) -> str:
    """Validate that file exists and has fasta extension."""
    if not value:
        return value
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"File not found: '{value}'")
    valid_ext = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz")
    if not value.lower().endswith(valid_ext):
        raise argparse.ArgumentTypeError(
            f"Expected FASTA file (.fa, .fasta), got: '{value}'"
        )
    return value


def existing_gtf_file(value) -> str:
    """Validate that file exists and has GTF extension."""
    if not value:
        return value
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"File not found: '{value}'")
    valid_ext = (".gtf", ".gtf.gz", ".gff", ".gff3")
    if not value.lower().endswith(valid_ext):
        raise argparse.ArgumentTypeError(
            f"Expected GTF/GFF file (.gtf, .gff), got: '{value}'"
        )
    return value


def existing_count_table(value) -> str:
    """Validate that file exists and has valid count table extension."""
    if not value:
        return value
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"File not found: '{value}'")
    valid_ext = (".csv", ".tsv", ".txt")
    if not value.lower().endswith(valid_ext):
        raise argparse.ArgumentTypeError(
            f"Expected count table (.csv, .tsv, .txt), got: '{value}'"
        )
    return value


def _range_errors(args: argparse.Namespace) -> list[str]:
    """Check min/max pairs and the k-mer length against the probe length."""
    errors = []

    if args.min_length > args.max_length:
        errors.append(
            f"--min-length ({args.min_length}) must be <= --max-length ({args.max_length})"
        )
    if args.min_tm > args.max_tm:
        errors.append(f"--min-tm ({args.min_tm}) must be <= --max-tm ({args.max_tm})")
    if args.min_gc > args.max_gc:
        errors.append(f"--min-gc ({args.min_gc}) must be <= --max-gc ({args.max_gc})")

    # Check kmer_length vs probe length
    if args.kmer_length >= args.min_length:
        errors.append(
            f"--kmer-length ({args.kmer_length}) must be < --min-length ({args.min_length})"
        )

    return errors


def _sequence_input_errors(args: argparse.Namespace) -> list[str]:
    """Check that a gene sequence source is given and NCBI names come in pairs."""
    errors = []

    # Check sequence input is provided (unless build_indices or analyze_probeset)
    if not args.build_indices and not args.analyze_probeset:
        has_sequence = (
            args.sequence_file
            or args.ensembl_id
            or (args.gene_name and args.organism_name)
        )
        if not has_sequence:
            errors.append(
                "Must provide gene sequence via one of:\n"
                "  --sequence-file FILE, or\n"
                "  --ensembl-id ID, or\n"
                "  --gene-name NAME --organism-name NAME"
            )

    # Check that gene_name requires organism_name and vice versa (for NCBI download)
    if args.sequence_file or args.ensembl_id:
        return errors
    if args.gene_name and not args.organism_name:
        errors.append("--gene-name requires --organism-name for NCBI download")
    if args.organism_name and not args.gene_name:
        errors.append("--organism-name requires --gene-name for NCBI download")

    return errors


# Options that only work together with --reference-annotation
_ANNOTATION_DEPENDENT_OPTIONS = {
    "encode_count_table": "--encode-count-table",
    "intergenic_off_targets": "--intergenic-off-targets",
    "filter_rrna": "--filter-rrna",
}


def _annotation_errors(args: argparse.Namespace) -> list[str]:
    """Check that options needing --reference-annotation have one."""
    return [
        f"{flag} requires --reference-annotation"
        for name, flag in _ANNOTATION_DEPENDENT_OPTIONS.items()
        if getattr(args, name, False) and not args.reference_annotation
    ]


def _lacks_exogenous_transcriptome(args: argparse.Namespace) -> bool:
    """Whether an exogenous design runs without transcriptome BLAST protection.

    Exogenous genes without transcriptome BLAST are under-protected:
    k-mer filtering is skipped, and genome alignment is fragile for short probes.
    Require explicit override to proceed without transcriptome.
    """
    if getattr(args, "is_endogenous", True):
        return False
    if args.reference_transcriptome or getattr(args, "allow_no_transcriptome", False):
        return False
    return not args.build_indices and not args.analyze_probeset


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate cross-parameter constraints after parsing."""
    errors = _range_errors(args) + _sequence_input_errors(args) + _annotation_errors(args)

    if _lacks_exogenous_transcriptome(args):
        errors.append(
            "Exogenous probe design requires --reference-transcriptome for reliable "
            "off-target detection. K-mer filtering is skipped for exogenous genes, so "
            "transcriptome BLAST is the primary off-target check.\n"
            "  To proceed without a transcriptome (reduced protection), add "
            "--allow-no-transcriptome."
        )

    if errors:
        parser.error("\n  " + "\n  ".join(errors))
