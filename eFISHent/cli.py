"""Command line interface."""

import argparse
import configparser
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

# Silence some troublemakers before heavy imports
warnings.filterwarnings("ignore")
logging.getLogger("luigi-interface").setLevel(level=logging.CRITICAL)
logging.getLogger("matplotlib.font_manager").setLevel(level=logging.CRITICAL)

import luigi

from . import __version__
from .constants import CLI_SHORTFORM, CONFIG_CLASSES
from .util import UniCode

if TYPE_CHECKING:
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtie2Index, BuildBowtieIndex
    from .kmers import BuildJellyfishIndex


GROUP_DESCRIPTIONS = {
    "General": "General configuration that will be used for all tasks.",
    "Run": "Options that change the behavior of the workflow.",
    "Sequence": "Details about the sequences the probe design will be performed on.",
    "Probe": "Probe filtering and design options.",
}
REQUIRED_PARAMS = ["reference_genome"]


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


# Parameter-specific type validators
PARAM_VALIDATORS = {
    # File validators
    "reference_genome": existing_fasta_file,
    "reference_annotation": existing_gtf_file,
    "reference_transcriptome": existing_fasta_file,
    "sequence_file": existing_fasta_file,
    "encode_count_table": existing_count_table,
    "analyze_probeset": existing_fasta_file,
    # Positive integers
    "threads": positive_int,
    "min_length": positive_int,
    "max_length": positive_int,
    "kmer_length": positive_int,
    "max_kmers": positive_int,
    "max_homopolymer_length": non_negative_int,
    "optimization_time_limit": positive_int,
    # Non-negative integers
    "spacing": non_negative_int,
    "max_off_targets": non_negative_int,
    "max_transcriptome_off_targets": non_negative_int,
    "sequence_similarity": non_negative_int,
    # Percentages
    "min_gc": percentage,
    "max_gc": percentage,
    "formamide_concentration": percentage,
    "max_expression_percentage": percentage,
    "blast_identity_threshold": percentage,
    "off_target_min_tm": float,
    "min_blast_match_length": non_negative_int,
    "max_probes_per_off_target": non_negative_int,
    "custom_rdna_fasta": existing_fasta_file,
    "max_cpg_fraction": float,
}

# Metavar hints for cleaner help output
PARAM_METAVAR = {
    "reference_genome": "FILE",
    "reference_annotation": "FILE",
    "sequence_file": "FILE",
    "encode_count_table": "FILE",
    "analyze_probeset": "FILE",
    "output_dir": "DIR",
    "threads": "N",
    "min_length": "N",
    "max_length": "N",
    "spacing": "N",
    "min_tm": "°C",
    "max_tm": "°C",
    "min_gc": "%",
    "max_gc": "%",
    "formamide_concentration": "%",
    "na_concentration": "mM",
    "max_off_targets": "N",
    "max_expression_percentage": "%",
    "kmer_length": "N",
    "max_kmers": "N",
    "max_deltag": "kcal/mol",
    "sequence_similarity": "%",
    "optimization_time_limit": "SEC",
    "gene_name": "NAME",
    "organism_name": "NAME",
    "ensembl_id": "ID",
    "optimization_method": "METHOD",
    "aligner": "ALIGNER",
    "max_homopolymer_length": "N",
    "max_transcriptome_off_targets": "N",
    "blast_identity_threshold": "%",
    "reference_transcriptome": "FILE",
    "off_target_min_tm": "°C",
    "min_blast_match_length": "N",
    "max_probes_per_off_target": "N",
    "custom_rdna_fasta": "FILE",
    "max_cpg_fraction": "FRAC",
}
_BOOL_METAVAR = "yes/no"
for _bool_param in [
    "filter_low_complexity", "build_indices", "save_intermediates",
    "is_plus_strand", "is_endogenous", "no_alternative_loci",
    "mask_repeats", "intergenic_off_targets", "filter_rrna",
    "adaptive_length", "filter_rdna_45s", "reject_cross_hybridization",
    "allow_no_transcriptome", "accessibility_scoring",
]:
    PARAM_METAVAR[_bool_param] = _BOOL_METAVAR

# Fixed choices for specific parameters
PARAM_CHOICES = {
    "optimization_method": ["greedy", "optimal"],
    "aligner": ["bowtie", "bowtie2"],
}


def get_parameter_type(param: luigi.Parameter) -> Any:
    """Get the type of a parameter."""
    if isinstance(param, luigi.IntParameter):
        return int
    elif isinstance(param, luigi.FloatParameter):
        return float
    elif isinstance(param, luigi.BoolParameter):
        return string_to_bool
    return str


def _add_utilities(parser: argparse.ArgumentParser) -> None:
    """Add the utility arguments to the parser."""
    utility = parser.add_argument_group("General utilities")
    utility.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this message.",
    )
    utility.add_argument(
        "-V",
        "--version",
        action="version",
        version="%(prog)s " + str(__version__),
        help="Show %(prog)s's version number.",
    )
    utility.add_argument(
        "-s",
        "--silent",
        action="store_true",
        help="Change the program output to silent hiding information on progress.",
    )
    utility.add_argument(
        "--check",
        action="store_true",
        help="Check if all external dependencies are installed and exit.",
    )
    utility.add_argument(
        "--update",
        action="store_true",
        help="Update eFISHent to the latest version.",
    )
    utility.add_argument(
        "--preset",
        type=str,
        default=None,
        metavar="NAME",
        help=(
            "Apply a parameter preset for common FISH protocols. "
            "Explicit arguments override preset values. "
            "Use --preset list to see available presets. "
            "[options: smfish, merfish, dna-fish, strict, relaxed, exogenous]"
        ),
    )
    utility.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    genome_group = parser.add_argument_group("Pre-built genome indices")
    genome_group.add_argument(
        "--genome",
        type=str,
        default=None,
        metavar="NAME",
        help=(
            "Use pre-built genome indices from Hugging Face Hub. "
            "Skips all index building. Mutually exclusive with --reference-genome. "
            "[options: hg38, mm39, dm6, human, mouse, fly]"
        ),
    )
    genome_group.add_argument(
        "--download-genome",
        type=str,
        default=None,
        metavar="NAME",
        help="Pre-download genome indices without running the pipeline.",
    )
    genome_group.add_argument(
        "--list-genomes",
        action="store_true",
        help="List available pre-built genomes and exit.",
    )
    genome_group.add_argument(
        "--index-cache-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="Override default cache directory for genome indices (~/.local/efishent/indices/).",
    )
    genome_group.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download genome indices even if already cached.",
    )


def _format_default(param: luigi.Parameter) -> str:
    """Format a parameter's default value for the help text."""
    if param._default is None or param._default == "":
        return ""
    if isinstance(param, luigi.BoolParameter):
        return f" (default: {'yes' if param._default else 'no'})"
    return f" (default: {param._default})"


def _argument_kwargs(name: str, param: luigi.Parameter) -> dict[str, Any]:
    """Build the add_argument keyword arguments for a luigi parameter."""
    # Use custom validator if available, otherwise default type
    if name in PARAM_VALIDATORS:
        param_type = PARAM_VALIDATORS[name]
    else:
        param_type = get_parameter_type(param)

    kwargs = {
        "type": param_type,
        "required": name in REQUIRED_PARAMS,
        "default": param._default,
        "help": f"{param.description}{_format_default(param)}",
    }

    # Add metavar for cleaner help output
    if name in PARAM_METAVAR:
        kwargs["metavar"] = PARAM_METAVAR[name]

    if name in PARAM_CHOICES:
        kwargs["choices"] = PARAM_CHOICES[name]
        kwargs["type"] = str  # Override type when using choices

    return kwargs


def _add_group(group: argparse._ArgumentGroup, config_class: luigi.Config) -> None:
    """Add a single configuration class/group to a parser."""
    for name, param in config_class().get_params():  # type: ignore
        kwargs = _argument_kwargs(name, param)
        group.add_argument(
            f"-{CLI_SHORTFORM.get(name)}",
            f"--{name.replace('_', '-')}",
            **kwargs,
        )


def _add_groups(parser: argparse.ArgumentParser) -> None:
    """Add the main option groups to the parser."""
    groups = [
        parser.add_argument_group(f"{name} options", description=description)
        for name, description in GROUP_DESCRIPTIONS.items()
    ]

    for group, config_class in zip(groups, CONFIG_CLASSES):
        _add_group(group, config_class)


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


def _efishent_home() -> Path:
    """Return eFISHent's per-user install directory (~/.local/efishent)."""
    return Path.home() / ".local" / "efishent"


def _ensure_deps_on_path() -> None:
    """Auto-discover eFISHent deps installed by install.sh.

    If deps were installed to ~/.local/efishent/deps/bin (the default),
    add them to PATH so users don't need to source activate.sh manually.
    """
    default_deps = _efishent_home() / "deps"
    deps_bin = default_deps / "bin"
    deps_edirect = default_deps / "edirect"

    if deps_bin.is_dir() and str(deps_bin) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{deps_bin}:{os.environ.get('PATH', '')}"

    if deps_edirect.is_dir() and str(deps_edirect) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = f"{deps_edirect}:{os.environ.get('PATH', '')}"

    # Also set library path for jellyfish/glpk shared libs
    deps_lib = default_deps / "lib"
    if deps_lib.is_dir():
        if sys.platform == "darwin":
            lib_var = "DYLD_LIBRARY_PATH"
        else:
            lib_var = "LD_LIBRARY_PATH"
        if str(deps_lib) not in os.environ.get(lib_var, ""):
            os.environ[lib_var] = f"{deps_lib}:{os.environ.get(lib_var, '')}"


def _get_tool_version(_name: str, args: list, pattern: str = r"[0-9]+\.[0-9]+\.?[0-9]*") -> str:
    """Try to get a tool's version string."""
    import re

    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=5)
        output = result.stdout + result.stderr
        match = re.search(pattern, output)
        return match.group(0) if match else "installed"
    except Exception:
        return ""


# Dependency definitions: {name: (check_args, version_args, needed_for)}
DEPENDENCIES = {
    "bowtie": {
        "version_args": ["bowtie", "--version"],
        "needed_for": "probe alignment (legacy aligner)",
    },
    "bowtie2": {
        "version_args": ["bowtie2", "--version"],
        "needed_for": "probe alignment (default aligner)",
    },
    "jellyfish": {
        "version_args": ["jellyfish", "--version"],
        "needed_for": "k-mer counting",
    },
    "blastn": {
        "version_args": ["blastn", "-version"],
        "needed_for": "transcriptome off-target filtering (optional)",
    },
    "makeblastdb": {
        "version_args": ["makeblastdb", "-version"],
        "needed_for": "building BLAST databases for transcriptome filtering (optional)",
    },
    "dustmasker": {
        "version_args": ["dustmasker", "-version-full"],
        "needed_for": "repeat masking for off-target filtering (optional)",
    },
    "glpsol": {
        "version_args": ["glpsol", "--version"],
        "needed_for": "optimal optimization (optional if using greedy)",
    },
    "esearch": {
        "version_args": ["esearch", "-version"],
        "needed_for": "NCBI gene download (optional if providing sequence file)",
    },
    "gffread": {
        "version_args": ["gffread", "--version"],
        "needed_for": "building transcriptome from genome + GTF (optional)",
    },
}


_FOLD_PURPOSE = "secondary structure prediction"


def _dependency_status(name: str, info: dict) -> dict:
    """Look up a single external tool and report where it is found."""
    path = shutil.which(name)
    if not path:
        return {"found": False, "needed_for": info["needed_for"]}
    version = _get_tool_version(name, info["version_args"])
    return {
        "found": True,
        "path": path,
        "version": version or path,
        "needed_for": info["needed_for"],
    }


def _system_fold() -> str | None:
    """Return the RNAstructure Fold binary on PATH, if any."""
    from .secondary_structure import _is_rnastructure_fold

    system_fold = shutil.which("Fold")
    if system_fold and _is_rnastructure_fold(system_fold):
        return system_fold
    return None


def _bundled_fold() -> Path | None:
    """Return the executable Fold binary bundled for this platform, if any."""
    fold_path = Path(__file__).resolve().parent
    if sys.platform.startswith("linux"):
        fold_bin = fold_path / "Fold_linux"
    elif sys.platform == "darwin":
        fold_bin = fold_path / "Fold_osx"
    else:
        return None

    if fold_bin.exists() and os.access(fold_bin, os.X_OK):
        return fold_bin
    return None


def _fold_status() -> dict:
    """Check for Fold binary: prefer system install, fall back to bundled."""
    system_fold = _system_fold()
    if system_fold:
        return {
            "found": True,
            "version": "system",
            "path": system_fold,
            "needed_for": _FOLD_PURPOSE,
        }

    fold_bin = _bundled_fold()
    if fold_bin:
        return {
            "found": True,
            "version": "bundled",
            "path": str(fold_bin),
            "needed_for": _FOLD_PURPOSE,
        }

    return {"found": False, "needed_for": _FOLD_PURPOSE}


def check_all_dependencies() -> dict[str, dict]:
    """Check all external dependencies and return status dict."""
    results = {
        name: _dependency_status(name, info) for name, info in DEPENDENCIES.items()
    }
    results["Fold"] = _fold_status()
    return results


def _needs_ncbi_download(args: argparse.Namespace) -> bool:
    """Whether the gene sequence has to be downloaded from NCBI."""
    if args.build_indices or getattr(args, "sequence_file", None):
        return False
    return bool(getattr(args, "gene_name", None) or getattr(args, "ensembl_id", None))


def _required_tools(args: argparse.Namespace) -> list[str]:
    """List the external tools needed for the current run mode (excluding Fold)."""
    # Aligner depends on --aligner setting
    aligner = getattr(args, "aligner", "bowtie2")
    required = ["bowtie2" if aligner == "bowtie2" else "bowtie", "jellyfish"]

    # BLAST required if transcriptome filtering, repeat masking, or rDNA filter is enabled
    if getattr(args, "reference_transcriptome", "") or getattr(args, "filter_rdna_45s", True):
        required += ["blastn", "makeblastdb"]
    if getattr(args, "mask_repeats", False):
        required.append("dustmasker")

    # glpsol only needed for optimal optimization
    if getattr(args, "optimization_method", "greedy") == "optimal":
        required.append("glpsol")

    # esearch/efetch only needed for NCBI download
    if _needs_ncbi_download(args):
        required.append("esearch")

    return required


def check_required_dependencies(args: argparse.Namespace) -> list[str]:
    """Check only the dependencies required for the current run mode.

    Returns list of missing dependency names (empty = all good).
    """
    # Fold is required for the full pipeline (not for index building)
    if not args.build_indices and not (_system_fold() or _bundled_fold()):
        return ["Fold (not found in PATH and bundled binary missing)"]

    return [name for name in _required_tools(args) if not shutil.which(name)]


def _build_parser() -> argparse.ArgumentParser:
    """Create the argument parser with rich help formatting."""
    from rich_argparse import RichHelpFormatter

    # Configure RichHelpFormatter for better alignment
    RichHelpFormatter.styles["argparse.metavar"] = "cyan italic"
    RichHelpFormatter.styles["argparse.args"] = "green"
    RichHelpFormatter.styles["argparse.groups"] = "bold blue"
    RichHelpFormatter.group_name_formatter = str.upper

    # Get terminal width, default to 120 if not available
    terminal_width = shutil.get_terminal_size((120, 24)).columns

    parser = argparse.ArgumentParser(
        prog="eFISHent",
        description=f"eFISHent v{__version__} {UniCode.fishing} {UniCode.dna} - RNA FISH probe designer",
        epilog=(
            'See the wiki at https://github.com/BBQuercus/eFISHent/wiki for details.\n'
            f"We hope you enjoy using eFISHent {UniCode.party}!"
        ),
        add_help=False,
        formatter_class=lambda prog: RichHelpFormatter(
            prog, max_help_position=40, width=terminal_width
        ),
    )
    _add_groups(parser)
    _add_utilities(parser)
    return parser


def _argv_value(flag: str) -> str | None:
    """Return the token following ``flag`` in sys.argv, if there is one."""
    idx = sys.argv.index(flag)
    if idx + 1 < len(sys.argv):
        return sys.argv[idx + 1]
    return None


# Utility flags handled before full validation (they don't need other args)
_STANDALONE_FLAGS = {
    "--check": "check",
    "--update": "update",
    "--list-genomes": "list_genomes",
}


def _parse_standalone_args(
    parser: argparse.ArgumentParser,
) -> argparse.Namespace | None:
    """Parse utility flags that bypass full argument parsing, if present."""
    for flag, dest in _STANDALONE_FLAGS.items():
        if flag in sys.argv:
            return argparse.Namespace(**{dest: True})

    if "--download-genome" not in sys.argv:
        return None

    genome_name = _argv_value("--download-genome")
    if genome_name is None:
        parser.error("--download-genome requires a genome name")
    cache_dir = None
    if "--index-cache-dir" in sys.argv:
        cache_dir = _argv_value("--index-cache-dir")
    return argparse.Namespace(
        download_genome=genome_name,
        index_cache_dir=cache_dir,
        force_download="--force-download" in sys.argv,
    )


def _make_reference_genome_optional(parser: argparse.ArgumentParser) -> None:
    """Drop the -g requirement (--genome sets it automatically later in main())."""
    for action in parser._actions:
        if "--reference-genome" in getattr(action, "option_strings", []):
            action.required = False
            break


def _apply_preset(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Apply preset defaults — explicit CLI args take precedence."""
    from .presets import PRESETS

    if args.preset not in PRESETS:
        from .presets import get_preset_names
        parser.error(
            f"Unknown preset '{args.preset}'. "
            f"Available: {', '.join(get_preset_names())}"
        )
    preset_params = PRESETS[args.preset]["params"]
    # Only apply preset values for params the user didn't explicitly set
    # We detect this by checking if the value matches the parser default
    for param_name, preset_value in preset_params.items():
        parser_default = parser.get_default(param_name)
        current_value = getattr(args, param_name, None)
        if current_value == parser_default:
            setattr(args, param_name, preset_value)


def _parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = _build_parser()
    try:
        if len(sys.argv) == 1:
            parser.print_help()
            parser.exit(0)
    except Exception as e:
        print(e)

    standalone_args = _parse_standalone_args(parser)
    if standalone_args is not None:
        return standalone_args

    # Handle --preset list before full parsing
    if "--preset" in sys.argv and _argv_value("--preset") == "list":
        from .presets import format_preset_list
        print(format_preset_list())
        sys.exit(0)

    # When --genome, --list-genomes, or --download-genome is used, -g is not required
    genome_flags = {"--genome", "--list-genomes", "--download-genome"}
    if genome_flags & set(sys.argv):
        _make_reference_genome_optional(parser)

    args = parser.parse_args()
    if args.preset:
        _apply_preset(args, parser)

    validate_args(args, parser)
    return args


def create_custom_config(args: argparse.Namespace, config_file: str) -> None:
    """Create a custom config file."""
    config = configparser.ConfigParser()
    config_path = Path(__file__).resolve().parent.joinpath("luigi.cfg").as_posix()
    config.read(config_path)

    for section, config_class in zip(
        ["GeneralConfig", "RunConfig", "SequenceConfig", "ProbeConfig"], CONFIG_CLASSES
    ):
        for name in config_class().get_param_names():  # type: ignore
            value = vars(args).get(name)
            config.set(section, name, str(value))
            if name == "threads":
                threads = min(vars(args).get(name), os.cpu_count())  # type: ignore
                config.set(section, name, str(threads))

    with open(config_file, "w") as f:
        config.write(f)


def set_logging_level(silent: bool, debug: bool) -> logging.Logger:
    """Set the logging level of luigi and custom logger."""
    from .console import get_rich_handler, set_silent_mode

    set_silent_mode(silent)
    luigi_level = "WARNING"
    handlers: list = []

    if debug:
        luigi_level = "DEBUG"
        custom_level = logging.DEBUG
        # File handler for debug mode
        file_handler = logging.FileHandler("efishent.log")
        file_handler.setFormatter(
            logging.Formatter(
                "%(asctime)s %(levelname)-4s [%(name)s] "
                "%(filename)s %(funcName)s %(lineno)d / %(thread)d - %(message)s"
            )
        )
        handlers.append(file_handler)
        # Also use Rich handler for console
        handlers.append(get_rich_handler())
    elif silent:
        custom_level = logging.WARNING
        # Use standard handler in silent mode
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        handlers.append(stream_handler)
    else:
        custom_level = logging.INFO
        # Use Rich handler for nice console output
        handlers.append(get_rich_handler())

    logging.basicConfig(handlers=handlers, force=True)  # type: ignore
    logging.getLogger("luigi").setLevel(luigi_level)
    # Keep luigi-interface at CRITICAL to suppress the execution summary
    # (eFISHent already shows its own completion panel).
    logging.getLogger("luigi-interface").setLevel(logging.CRITICAL)
    luigi.interface.core.log_level = luigi_level  # type: ignore

    logger = logging.getLogger("custom-logger")
    logger.setLevel(custom_level)
    return logger


_PYPI_JSON_URL = "https://pypi.org/pypi/efishent/json"
_UPDATE_CHECK_INTERVAL = 86400  # 24 hours


def _notify_if_newer(latest: str) -> None:
    """Print an update notice if ``latest`` is newer than the running version."""
    if not latest or latest == __version__:
        return

    from packaging.version import Version

    if Version(latest) > Version(__version__):
        print(
            f"\n  Update available: v{__version__} -> v{latest}"
            f"  (run 'efishent --update' to upgrade)\n"
        )


def _cached_latest_version(cache_file: Path) -> str | None:
    """Return the cached latest version, or None if the cache is missing or stale."""
    if not cache_file.exists():
        return None
    with open(cache_file) as f:
        cache = json.load(f)
    last_check = cache.get("timestamp", 0)
    if time.time() - last_check >= _UPDATE_CHECK_INTERVAL:
        return None
    return cache.get("latest_version", "")


def _https_get_json(url: str) -> dict:
    """Fetch a JSON document over https; other URL schemes are refused."""
    import urllib.parse
    import urllib.request

    if urllib.parse.urlparse(url).scheme != "https":
        raise ValueError(f"Refusing to open non-https URL: {url}")

    # Only https handlers are installed, so redirects to other schemes fail too
    opener = urllib.request.OpenerDirector()
    for handler in (
        urllib.request.ProxyHandler(),
        urllib.request.HTTPSHandler(),
        urllib.request.HTTPRedirectHandler(),
        urllib.request.HTTPDefaultErrorHandler(),
        urllib.request.HTTPErrorProcessor(),
        urllib.request.UnknownHandler(),
    ):
        opener.add_handler(handler)

    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with opener.open(req, timeout=5) as resp:
        return json.loads(resp.read().decode())


def _run_update_check(cache_file: Path) -> None:
    """Notify about a newer release, querying PyPI unless a fresh cache exists."""
    try:
        # Check cache first (skip if checked within last 24 hours)
        latest = _cached_latest_version(cache_file)
        if latest is None:
            data = _https_get_json(_PYPI_JSON_URL)
            latest = data.get("info", {}).get("version", "")

            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w") as f:
                json.dump({"timestamp": time.time(), "latest_version": latest}, f)

        _notify_if_newer(latest)
    except Exception as e:
        # Never fail the pipeline for an update check
        logging.getLogger("custom-logger").debug(f"Update check failed: {e}")


def check_for_updates() -> None:
    """Check PyPI for a newer version and notify the user.

    Non-blocking: runs in a background thread so it doesn't slow startup.
    Caches the result for 24 hours to avoid repeated network requests.
    """
    import threading

    cache_file = _efishent_home() / "update_check.json"
    thread = threading.Thread(target=_run_update_check, args=(cache_file,), daemon=True)
    thread.start()


_INSTALL_SCRIPT_URL = (
    "https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh"
)


def _upgrade_package(python: str) -> None:
    """Upgrade the eFISHent package in place, exiting on failure."""
    if shutil.which("uv"):
        cmd = ["uv", "pip", "install", "--python", python, "--upgrade", "efishent"]
    else:
        cmd = [python, "-m", "pip", "install", "--upgrade", "efishent"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except Exception as e:
        print(f"Package update failed: {e}")
        sys.exit(1)
    if result.returncode != 0:
        print(f"Package update failed:\n{result.stderr.strip()}")
        sys.exit(1)


def _installed_version(python: str) -> str:
    """Return the eFISHent version importable by ``python``."""
    try:
        new_version_result = subprocess.run(
            [python, "-c", "from eFISHent import __version__; print(__version__)"],
            capture_output=True, text=True, timeout=10,
        )
        return new_version_result.stdout.strip()
    except Exception:
        return "unknown"


def _download_command(url: str) -> list[str] | None:
    """Return a curl/wget command that writes ``url`` to stdout, if available."""
    if shutil.which("curl"):
        return ["curl", "-fsSL", url]
    if shutil.which("wget"):
        return ["wget", "-qO-", url]
    return None


def _run_install_script(dl_cmd: list[str], prefix: Path) -> None:
    """Download install.sh and re-run it with --deps-only to update/verify deps."""
    dl_result = subprocess.run(dl_cmd, capture_output=True, text=True, timeout=30)
    if dl_result.returncode != 0:
        print("Warning: could not download install script.")
        return
    result = subprocess.run(
        ["sh", "-s", "--", "--deps-only", "--no-modify-rc",
         "--prefix", str(prefix)],
        input=dl_result.stdout, text=True, timeout=600,
    )
    if result.returncode != 0:
        print("Warning: dependency update had issues.")


def _update_dependencies(prefix: Path) -> None:
    """Update dependencies via install.sh if eFISHent was installed with it."""
    if not prefix.is_dir():
        print(f"No install directory at {prefix}, skipping dependency update.")
        print("Run the full installer to set up dependencies:")
        print(f"  curl -LsSf {_INSTALL_SCRIPT_URL} | sh")
        return

    dl_cmd = _download_command(_INSTALL_SCRIPT_URL)
    if not dl_cmd:
        print("Warning: neither curl nor wget found, skipping dependency update.")
        return

    try:
        _run_install_script(dl_cmd, prefix)
    except Exception as e:
        print(f"Warning: dependency update skipped ({e})")


def _print_dependency_summary() -> None:
    """Print which dependencies are found and how to install missing ones."""
    _ensure_deps_on_path()
    results = check_all_dependencies()
    all_found = all(r["found"] for r in results.values())
    for name, info in results.items():
        if info["found"]:
            version = info.get("version", "")
            print(f"  OK  {name}: {version}")
        else:
            print(f"  MISS {name}: not found ({info['needed_for']})")

    if not all_found:
        print("\nSome dependencies are missing. Run the full installer:")
        print(f"  curl -LsSf {_INSTALL_SCRIPT_URL} | sh")


def self_update() -> None:
    """Update eFISHent and its dependencies to the latest version."""
    python = sys.executable
    current_version = __version__

    print(f"Current version: {current_version}")

    # Step 1: Update the Python package
    print("\nUpdating eFISHent package...")
    _upgrade_package(python)

    new_version = _installed_version(python)
    if new_version == current_version:
        print(f"Already up to date (v{current_version}).")
    else:
        print(f"Updated: v{current_version} -> v{new_version}")

    # Step 2: Update dependencies via install.sh
    print("\nChecking dependencies...")
    _update_dependencies(_efishent_home())

    # Step 3: Verify dependencies
    print()
    _print_dependency_summary()


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{minutes}m {secs}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"


def _prefill_cached_steps(add_cached_step_fn) -> None:
    """Check which pipeline steps are already cached and add them to the display."""
    from .alignment import AlignProbeCandidates
    from .basic_filtering import BasicFiltering
    from .cleanup import CleanUpOutput
    from .generate_probes import GenerateAllProbes
    from .kmers import KMerFiltering
    from .optimization import OptimizeProbeCoverage

    # Map stage names to Luigi task classes
    from .prepare_sequence import PrepareSequence
    from .rdna_filter import FilterRibosomalRNA
    from .secondary_structure import SecondaryStructureFiltering
    from .transcriptome_filter import TranscriptomeFiltering
    from .util import _get_pipeline_stages

    task_classes = {
        "PrepareSequence": PrepareSequence,
        "GenerateAllProbes": GenerateAllProbes,
        "BasicFiltering": BasicFiltering,
        "FilterRibosomalRNA": FilterRibosomalRNA,
        "AlignProbeCandidates": AlignProbeCandidates,
        "TranscriptomeFiltering": TranscriptomeFiltering,
        "KMerFiltering": KMerFiltering,
        "SecondaryStructureFiltering": SecondaryStructureFiltering,
        "OptimizeProbeCoverage": OptimizeProbeCoverage,
        "CleanUpOutput": CleanUpOutput,
    }

    # Get config-aware stage definitions and sort by step number
    stages = _get_pipeline_stages()
    numbered_steps = sorted(
        [(v["order"], v["total"], v["desc"], task_classes.get(k))
         for k, v in stages.items() if v["order"] > 0],
        key=lambda x: x[0],
    )

    for order, total, desc, task_cls in numbered_steps:
        try:
            if task_cls().complete():
                add_cached_step_fn(order, total, desc)
            else:
                break  # Stop at first incomplete step
        except Exception:
            break


def _run_dependency_check() -> None:
    """Handle --check: show dependency report and exit."""
    from .console import print_dependency_check, print_header

    print_header(__version__)
    results = check_all_dependencies()
    print_dependency_check(results)
    sys.exit(0 if all(r["found"] for r in results.values()) else 1)


def _list_genomes() -> None:
    """Handle --list-genomes: show available genomes and exit."""
    from .prebuilt import list_available_genomes

    print("Available pre-built genomes:")
    for genome in list_available_genomes():
        print(f"  {genome['display']}")
    print()
    print("Usage: efishent --genome hg38 --gene-name TP53 --organism-name 'Homo sapiens'")
    print("Pre-download: efishent --download-genome hg38")
    sys.exit(0)


def _download_genome_only(args: argparse.Namespace, genome_name: str) -> None:
    """Handle --download-genome: pre-download and exit."""
    from .prebuilt import download_genome, resolve_genome

    genome_id = resolve_genome(genome_name)
    cache_dir = getattr(args, "index_cache_dir", None)
    force = getattr(args, "force_download", False)
    download_genome(genome_id, cache_dir=cache_dir, force=force)
    print(f"Genome {genome_id} downloaded successfully.")
    sys.exit(0)


def _run_standalone_modes(args: argparse.Namespace) -> None:
    """Run utility modes that exit without running the pipeline."""
    if getattr(args, "check", False):
        _run_dependency_check()

    # Handle --update: self-update and exit
    if getattr(args, "update", False):
        self_update()
        sys.exit(0)

    if getattr(args, "list_genomes", False):
        _list_genomes()

    download_genome_arg = getattr(args, "download_genome", None)
    if download_genome_arg:
        _download_genome_only(args, download_genome_arg)


def _fetch_supplemental_kmer_index(
    genome_id: str, kmer_length: int, cache_dir: str | None
) -> None:
    """Auto-fetch supplemental jellyfish index for non-default kmer_length.

    genome_15.jf is always bundled; other kl values are hosted separately
    on HuggingFace and downloaded here on first use.
    """
    from .console import console
    from .prebuilt import try_download_supplemental_index

    console.print(
        f"  Checking HuggingFace for pre-built {kmer_length}-mer index...",
        end="",
    )
    found = try_download_supplemental_index(genome_id, kmer_length, cache_dir)
    if found:
        console.print(" [green]found[/green]")
    else:
        console.print(" [yellow]not available — will build locally[/yellow]")


def _use_prebuilt_genome(args: argparse.Namespace, genome_name: str) -> None:
    """Handle --genome: resolve and set reference paths."""
    from .prebuilt import (
        download_genome,
        get_reference_paths,
        is_genome_cached,
        resolve_genome,
    )

    if getattr(args, "reference_genome", None):
        print("Error: --genome and --reference-genome are mutually exclusive.")
        sys.exit(1)

    genome_id = resolve_genome(genome_name)
    cache_dir = getattr(args, "index_cache_dir", None)

    if not is_genome_cached(genome_id, cache_dir):
        print(f"Genome {genome_id} not cached. Downloading...")
        download_genome(genome_id, cache_dir=cache_dir)

    if args.kmer_length != 15:
        _fetch_supplemental_kmer_index(genome_id, args.kmer_length, cache_dir)

    paths = get_reference_paths(genome_id, cache_dir)
    args.reference_genome = paths["reference_genome"]
    args.reference_annotation = paths["reference_annotation"]
    args.reference_transcriptome = paths["reference_transcriptome"]


def _announce_start(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Print header banner (and check for updates) or log the start when silent."""
    if args.silent:
        logger.info(f"{UniCode.fishing} eFISHent v{__version__} starting...")
        return

    from .console import print_header

    print_header(__version__)
    check_for_updates()


def _exit_if_missing_dependencies(
    args: argparse.Namespace, logger: logging.Logger
) -> None:
    """Check required dependencies before starting the pipeline."""
    missing = check_required_dependencies(args)
    if not missing:
        return

    if not args.silent:
        from .console import print_missing_deps_error

        print_missing_deps_error(missing)
    else:
        logger.error(
            f"Missing required dependencies: {', '.join(missing)}. "
            "Run efishent --check for details."
        )
    sys.exit(1)


def _report_parameter_warnings(
    args: argparse.Namespace, logger: logging.Logger
) -> None:
    """Show pre-flight parameter warnings."""
    param_warnings = validate_parameter_warnings(args)
    if not param_warnings:
        return

    if not args.silent:
        from .console import print_parameter_warnings

        print_parameter_warnings(param_warnings)
    else:
        for w in param_warnings:
            logger.warning(w.replace("\n", " "))


def _task_classes(args: argparse.Namespace) -> list[type[luigi.Task]]:
    """Select the luigi tasks to run for the current mode."""
    # Lazy imports - only load heavy modules after arg parsing
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtie2Index, BuildBowtieIndex
    from .kmers import BuildJellyfishIndex

    if args.build_indices:
        if args.aligner == "bowtie2":
            index_task = BuildBowtie2Index
        else:
            index_task = BuildBowtieIndex
        return [BuildJellyfishIndex, index_task]
    if args.analyze_probeset:
        return [AnalyzeProbeset]
    return [CleanUpOutput]


def _run_tasks(tasks: list[luigi.Task], args: argparse.Namespace) -> None:
    """Run the tasks, with progress tracking unless silent."""
    if args.silent:
        luigi.build(tasks, local_scheduler=True)
        return

    from .console import add_cached_step, pipeline_progress

    if args.analyze_probeset:
        # Analysis has 7 steps
        with pipeline_progress(total_stages=7, mode="analyze"):
            luigi.build(tasks, local_scheduler=True)
    elif args.build_indices:
        with pipeline_progress(total_stages=2, mode="pipeline"):
            luigi.build(tasks, local_scheduler=True)
    else:
        from .util import _get_pipeline_total
        with pipeline_progress(total_stages=_get_pipeline_total()):
            _prefill_cached_steps(add_cached_step)
            luigi.build(tasks, local_scheduler=True)


def _task_output_files(task: luigi.Task) -> list[str]:
    """Collect output file paths of a task."""
    task_output = task.output()
    if isinstance(task_output, dict):
        return [t.path for t in task_output.values() if hasattr(t, "path")]
    if hasattr(task_output, "path"):
        return [task_output.path]
    return []


def _report_completion(
    task: luigi.Task,
    args: argparse.Namespace,
    logger: logging.Logger,
    duration: float,
) -> None:
    """Show the completion summary for the final task."""
    output_files = _task_output_files(task)

    if args.silent:
        logger.info(
            f"{UniCode.party} eFISHent has finished running in {format_duration(duration)}!"
        )
        for f in output_files:
            logger.info(f"  Output: {f}")
        return

    from .console import print_completion

    # Get summary stats and probe data from CleanUpOutput task if available
    print_completion(
        format_duration(duration),
        output_files,
        getattr(task, "_summary", None),
        getattr(task, "_probe_df", None),
        getattr(task, "_verification", None),
    )


def main():
    """Run eFISHent tasks."""
    # Auto-discover deps installed by install.sh
    _ensure_deps_on_path()

    start_time = time.time()
    args = _parse_args()

    _run_standalone_modes(args)

    genome_arg = getattr(args, "genome", None)
    if genome_arg:
        _use_prebuilt_genome(args, genome_arg)

    # Validate that reference_genome is set (either via -g or --genome)
    if not getattr(args, "reference_genome", None):
        print("eFISHent: error: the following arguments are required: -g/--reference-genome (or use --genome)")
        sys.exit(2)

    logger = set_logging_level(args.silent, args.debug)
    _announce_start(args, logger)
    _exit_if_missing_dependencies(args, logger)
    _report_parameter_warnings(args, logger)

    task_classes = _task_classes(args)

    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        create_custom_config(args, config_file)
        luigi.configuration.add_config_path(config_file)

        tasks = [task_cls() for task_cls in task_classes]
        _run_tasks(tasks, args)

    duration = time.time() - start_time

    if tasks[-1].complete():
        _report_completion(tasks[-1], args, logger, duration)
