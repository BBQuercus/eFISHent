"""Command line interface."""

from pathlib import Path
from typing import Any, Dict, List, TYPE_CHECKING
import argparse
import configparser
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

# Silence some troublemakers before heavy imports
warnings.filterwarnings("ignore")
logging.getLogger("luigi-interface").setLevel(level=logging.CRITICAL)
logging.getLogger("matplotlib.font_manager").setLevel(level=logging.CRITICAL)

import luigi

from . import __version__
from .constants import CLI_SHORTFORM
from .constants import CONFIG_CLASSES
from .util import UniCode

if TYPE_CHECKING:
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtieIndex, BuildBowtie2Index
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
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'")
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"Value must be >= 1, got {ivalue}")
    return ivalue


def non_negative_int(value) -> int:
    """Validate that value is a non-negative integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'")
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"Value must be >= 0, got {ivalue}")
    return ivalue


def percentage(value) -> float:
    """Validate that value is a percentage (0-100)."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid number: '{value}'")
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


def _add_group(group: argparse._ArgumentGroup, config_class: luigi.Config) -> None:
    """Add a single configuration class/group to a parser."""
    for name, param in config_class().get_params():  # type: ignore
        # Use custom validator if available, otherwise default type
        if name in PARAM_VALIDATORS:
            param_type = PARAM_VALIDATORS[name]
        else:
            param_type = get_parameter_type(param)

        is_required = name in REQUIRED_PARAMS

        # Format default value display
        if param._default is None or param._default == "":
            default_str = ""
        elif isinstance(param, luigi.BoolParameter):
            default_str = f" (default: {'yes' if param._default else 'no'})"
        else:
            default_str = f" (default: {param._default})"

        # Build argument kwargs
        kwargs = {
            "type": param_type,
            "required": is_required,
            "default": param._default,
            "help": f"{param.description}{default_str}",
        }

        # Add metavar for cleaner help output
        if name in PARAM_METAVAR:
            kwargs["metavar"] = PARAM_METAVAR[name]

        # Add choices for specific parameters
        if name == "optimization_method":
            kwargs["choices"] = ["greedy", "optimal"]
            kwargs["type"] = str  # Override type when using choices

        if name == "aligner":
            kwargs["choices"] = ["bowtie", "bowtie2"]
            kwargs["type"] = str

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


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate cross-parameter constraints after parsing."""
    errors = []

    # Check min/max pairs
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
    if (
        args.gene_name
        and not args.organism_name
        and not args.sequence_file
        and not args.ensembl_id
    ):
        errors.append("--gene-name requires --organism-name for NCBI download")
    if (
        args.organism_name
        and not args.gene_name
        and not args.sequence_file
        and not args.ensembl_id
    ):
        errors.append("--organism-name requires --gene-name for NCBI download")

    # Check encode_count_table requires reference_annotation
    if args.encode_count_table and not args.reference_annotation:
        errors.append("--encode-count-table requires --reference-annotation")

    # Check intergenic_off_targets requires reference_annotation
    if getattr(args, "intergenic_off_targets", False) and not args.reference_annotation:
        errors.append("--intergenic-off-targets requires --reference-annotation")

    # Check filter_rrna requires reference_annotation
    if getattr(args, "filter_rrna", False) and not args.reference_annotation:
        errors.append("--filter-rrna requires --reference-annotation")

    # Exogenous genes without transcriptome BLAST are under-protected:
    # k-mer filtering is skipped, and genome alignment is fragile for short probes.
    # Require explicit override to proceed without transcriptome.
    if (
        not getattr(args, "is_endogenous", True)
        and not args.reference_transcriptome
        and not getattr(args, "allow_no_transcriptome", False)
        and not args.build_indices
        and not args.analyze_probeset
    ):
        errors.append(
            "Exogenous probe design requires --reference-transcriptome for reliable "
            "off-target detection. K-mer filtering is skipped for exogenous genes, so "
            "transcriptome BLAST is the primary off-target check.\n"
            "  To proceed without a transcriptome (reduced protection), add "
            "--allow-no-transcriptome."
        )

    if errors:
        parser.error("\n  " + "\n  ".join(errors))


def validate_parameter_warnings(args: argparse.Namespace) -> List[str]:
    """Check for risky parameter combinations and return warning messages.

    These are non-fatal warnings — the pipeline will still run.
    """
    # Skip for non-pipeline modes
    if args.build_indices or args.analyze_probeset:
        return []

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

    # High-GC target warning: compute gene GC from sequence file if available
    seq_file = getattr(args, "sequence_file", None)
    if seq_file and os.path.isfile(seq_file):
        try:
            import Bio.SeqIO
            import Bio.SeqUtils
            seq = next(Bio.SeqIO.parse(seq_file, "fasta"))
            gene_gc = float(Bio.SeqUtils.gc_fraction(seq.seq)) * 100
            if gene_gc > 55:
                warnings.append(
                    f"High-GC target detected ({gene_gc:.1f}% GC) \u2014 elevated risk of "
                    "nuclear/nucleolar background.\n"
                    "  Probes from high-GC genes cross-hybridize with rRNA (~60-67% GC).\n"
                    "  Consider tightening --max-gc (e.g., 60) or using --max-cpg-fraction 0.10."
                )
        except Exception:
            pass

    return warnings


def _ensure_deps_on_path() -> None:
    """Auto-discover eFISHent deps installed by install.sh.

    If deps were installed to ~/.local/efishent/deps/bin (the default),
    add them to PATH so users don't need to source activate.sh manually.
    """
    default_deps = Path.home() / ".local" / "efishent" / "deps"
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


def _get_tool_version(name: str, args: list, pattern: str = r"[0-9]+\.[0-9]+\.?[0-9]*") -> str:
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


def check_all_dependencies() -> Dict[str, dict]:
    """Check all external dependencies and return status dict."""
    results = {}

    for name, info in DEPENDENCIES.items():
        path = shutil.which(name)
        if path:
            version = _get_tool_version(name, info["version_args"])
            results[name] = {
                "found": True,
                "path": path,
                "version": version or path,
                "needed_for": info["needed_for"],
            }
        else:
            results[name] = {
                "found": False,
                "needed_for": info["needed_for"],
            }

    # Check bundled Fold binary
    fold_path = Path(__file__).resolve().parent
    if sys.platform.startswith("linux"):
        fold_bin = fold_path / "Fold_linux"
    elif sys.platform == "darwin":
        fold_bin = fold_path / "Fold_osx"
    else:
        fold_bin = None

    if fold_bin and fold_bin.exists() and os.access(fold_bin, os.X_OK):
        results["Fold"] = {
            "found": True,
            "version": "bundled",
            "path": str(fold_bin),
            "needed_for": "secondary structure prediction",
        }
    else:
        results["Fold"] = {
            "found": False,
            "needed_for": "secondary structure prediction",
        }

    return results


def check_required_dependencies(args: argparse.Namespace) -> List[str]:
    """Check only the dependencies required for the current run mode.

    Returns list of missing dependency names (empty = all good).
    """
    # Aligner depends on --aligner setting
    aligner = getattr(args, "aligner", "bowtie2")
    if aligner == "bowtie2":
        required = ["bowtie2", "jellyfish"]
    else:
        required = ["bowtie", "jellyfish"]

    # BLAST required if transcriptome filtering, repeat masking, or rDNA filter is enabled
    if getattr(args, "reference_transcriptome", ""):
        required.append("blastn")
        required.append("makeblastdb")
    if getattr(args, "filter_rdna_45s", True) and "blastn" not in required:
        required.append("blastn")
        required.append("makeblastdb")
    if getattr(args, "mask_repeats", False):
        required.append("dustmasker")

    # Fold is required for the full pipeline (not for index building)
    if not args.build_indices:
        fold_path = Path(__file__).resolve().parent
        if sys.platform.startswith("linux"):
            fold_bin = fold_path / "Fold_linux"
        elif sys.platform == "darwin":
            fold_bin = fold_path / "Fold_osx"
        else:
            fold_bin = None

        if not (fold_bin and fold_bin.exists() and os.access(fold_bin, os.X_OK)):
            return ["Fold (bundled binary missing or not executable)"]

    # glpsol only needed for optimal optimization
    if getattr(args, "optimization_method", "greedy") == "optimal":
        required.append("glpsol")

    # esearch/efetch only needed for NCBI download
    needs_ncbi = (
        not args.build_indices
        and not getattr(args, "sequence_file", None)
        and (getattr(args, "gene_name", None) or getattr(args, "ensembl_id", None))
    )
    if needs_ncbi:
        required.append("esearch")

    missing = [name for name in required if not shutil.which(name)]
    return missing


def _parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
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
    try:
        if len(sys.argv) == 1:
            parser.print_help()
            parser.exit(0)
    except Exception as e:
        print(e)

    # Handle --check/--update before full validation (doesn't need other args)
    if "--check" in sys.argv:
        args = argparse.Namespace(check=True)
        return args

    if "--update" in sys.argv:
        args = argparse.Namespace(update=True)
        return args

    # Handle --preset list before full parsing
    if "--preset" in sys.argv:
        idx = sys.argv.index("--preset")
        if idx + 1 < len(sys.argv) and sys.argv[idx + 1] == "list":
            from .presets import format_preset_list
            print(format_preset_list())
            sys.exit(0)

    args = parser.parse_args()

    # Apply preset defaults — explicit CLI args take precedence
    if args.preset:
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
    logging.getLogger("luigi-interface").setLevel(luigi_level)
    luigi.interface.core.log_level = luigi_level  # type: ignore

    logger = logging.getLogger("custom-logger")
    logger.setLevel(custom_level)
    return logger


def self_update() -> None:
    """Update eFISHent and its dependencies to the latest version."""
    python = sys.executable
    current_version = __version__
    install_script_url = (
        "https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install.sh"
    )

    print(f"Current version: {current_version}")

    # Step 1: Update the Python package
    print("\nUpdating eFISHent package...")
    if shutil.which("uv"):
        cmd = ["uv", "pip", "install", "--python", python, "--upgrade", "efishent"]
    else:
        cmd = [python, "-m", "pip", "install", "--upgrade", "efishent"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"Package update failed:\n{result.stderr.strip()}")
            sys.exit(1)
    except Exception as e:
        print(f"Package update failed: {e}")
        sys.exit(1)

    # Check new version
    try:
        new_version_result = subprocess.run(
            [python, "-c", "from eFISHent import __version__; print(__version__)"],
            capture_output=True, text=True, timeout=10,
        )
        new_version = new_version_result.stdout.strip()
    except Exception:
        new_version = "unknown"

    if new_version == current_version:
        print(f"Already up to date (v{current_version}).")
    else:
        print(f"Updated: v{current_version} -> v{new_version}")

    # Step 2: Update dependencies via install.sh
    print("\nChecking dependencies...")
    default_prefix = Path.home() / ".local" / "efishent"

    if default_prefix.is_dir():
        # Re-run install.sh with --deps-only to update/verify deps
        dl_cmd = None
        if shutil.which("curl"):
            dl_cmd = ["curl", "-fsSL", install_script_url]
        elif shutil.which("wget"):
            dl_cmd = ["wget", "-qO-", install_script_url]

        if dl_cmd:
            try:
                dl_result = subprocess.run(dl_cmd, capture_output=True, text=True, timeout=30)
                if dl_result.returncode == 0:
                    result = subprocess.run(
                        ["sh", "-s", "--", "--deps-only", "--no-modify-rc",
                         "--prefix", str(default_prefix)],
                        input=dl_result.stdout, text=True, timeout=600,
                    )
                    if result.returncode != 0:
                        print("Warning: dependency update had issues.")
                else:
                    print("Warning: could not download install script.")
            except Exception as e:
                print(f"Warning: dependency update skipped ({e})")
        else:
            print("Warning: neither curl nor wget found, skipping dependency update.")
    else:
        print(f"No install directory at {default_prefix}, skipping dependency update.")
        print("Run the full installer to set up dependencies:")
        print(f"  curl -LsSf {install_script_url} | sh")

    # Step 3: Verify dependencies
    print()
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
        print(f"  curl -LsSf {install_script_url} | sh")


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


def main():
    """Run eFISHent tasks."""
    from .console import (
        pipeline_progress,
        print_completion,
        print_dependency_check,
        print_header,
        print_missing_deps_error,
        print_parameter_warnings,
    )

    # Auto-discover deps installed by install.sh
    _ensure_deps_on_path()

    start_time = time.time()
    args = _parse_args()

    # Handle --check: show dependency report and exit
    if getattr(args, "check", False):
        print_header(__version__)
        results = check_all_dependencies()
        print_dependency_check(results)
        sys.exit(0 if all(r["found"] for r in results.values()) else 1)

    # Handle --update: self-update and exit
    if getattr(args, "update", False):
        self_update()
        sys.exit(0)

    # Handle --list-genomes: show available genomes and exit
    if getattr(args, "list_genomes", False):
        from .prebuilt import list_available_genomes
        print("Available pre-built genomes:")
        for genome in list_available_genomes():
            print(f"  {genome['display']}")
        print()
        print("Usage: efishent --genome hg38 --gene-name TP53 --organism-name 'Homo sapiens'")
        print("Pre-download: efishent --download-genome hg38")
        sys.exit(0)

    # Handle --download-genome: pre-download and exit
    download_genome_arg = getattr(args, "download_genome", None)
    if download_genome_arg:
        from .prebuilt import download_genome, resolve_genome
        genome_id = resolve_genome(download_genome_arg)
        cache_dir = getattr(args, "index_cache_dir", None)
        force = getattr(args, "force_download", False)
        download_genome(genome_id, cache_dir=cache_dir, force=force)
        print(f"Genome {genome_id} downloaded successfully.")
        sys.exit(0)

    # Handle --genome: resolve and set reference paths
    genome_arg = getattr(args, "genome", None)
    if genome_arg:
        from .prebuilt import download_genome, get_reference_paths, is_genome_cached, resolve_genome
        if getattr(args, "reference_genome", None):
            print("Error: --genome and --reference-genome are mutually exclusive.")
            sys.exit(1)

        genome_id = resolve_genome(genome_arg)
        cache_dir = getattr(args, "index_cache_dir", None)

        if not is_genome_cached(genome_id, cache_dir):
            print(f"Genome {genome_id} not cached. Downloading...")
            download_genome(genome_id, cache_dir=cache_dir)

        paths = get_reference_paths(genome_id, cache_dir)
        args.reference_genome = paths["reference_genome"]
        args.reference_annotation = paths["reference_annotation"]
        args.reference_transcriptome = paths["reference_transcriptome"]

    logger = set_logging_level(args.silent, args.debug)

    # Print header banner
    if not args.silent:
        print_header(__version__)
    else:
        logger.info(f"{UniCode.fishing} eFISHent v{__version__} starting...")

    # Check required dependencies before starting the pipeline
    missing = check_required_dependencies(args)
    if missing:
        if not args.silent:
            print_missing_deps_error(missing)
        else:
            logger.error(
                f"Missing required dependencies: {', '.join(missing)}. "
                "Run efishent --check for details."
            )
        sys.exit(1)

    # Pre-flight parameter warnings
    param_warnings = validate_parameter_warnings(args)
    if param_warnings:
        if not args.silent:
            print_parameter_warnings(param_warnings)
        else:
            for w in param_warnings:
                logger.warning(w.replace("\n", " "))

    # Lazy imports - only load heavy modules after arg parsing
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtieIndex, BuildBowtie2Index
    from .kmers import BuildJellyfishIndex

    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        create_custom_config(args, config_file)
        luigi.configuration.add_config_path(config_file)

        tasks: List[luigi.Task] = []
        if args.build_indices:
            if args.aligner == "bowtie2":
                index_task = BuildBowtie2Index()
            else:
                index_task = BuildBowtieIndex()
            tasks = [BuildJellyfishIndex(), index_task]
        elif args.analyze_probeset:
            tasks = [AnalyzeProbeset()]
        else:
            tasks = [CleanUpOutput()]

        # Run with progress tracking
        if not args.silent:
            if args.analyze_probeset:
                # Analysis has 9 steps (preparing + 8 analysis plots)
                with pipeline_progress(total_stages=9, mode="analyze"):
                    luigi.build(tasks, local_scheduler=True)
            elif args.build_indices:
                with pipeline_progress(total_stages=2, mode="pipeline"):
                    luigi.build(tasks, local_scheduler=True)
            else:
                with pipeline_progress(total_stages=9):
                    luigi.build(tasks, local_scheduler=True)
        else:
            luigi.build(tasks, local_scheduler=True)

    duration = time.time() - start_time

    if tasks[-1].complete():
        # Collect output file paths
        output_files = []
        task_output = tasks[-1].output()
        if isinstance(task_output, dict):
            output_files = [t.path for t in task_output.values() if hasattr(t, "path")]
        elif hasattr(task_output, "path"):
            output_files = [task_output.path]

        # Get summary stats and probe data from CleanUpOutput task if available
        summary = getattr(tasks[-1], "_summary", None)
        probe_df = getattr(tasks[-1], "_probe_df", None)
        verification = getattr(tasks[-1], "_verification", None)

        if not args.silent:
            print_completion(
                format_duration(duration), output_files, summary, probe_df,
                verification,
            )
        else:
            logger.info(
                f"{UniCode.party} eFISHent has finished running in {format_duration(duration)}!"
            )
            for f in output_files:
                logger.info(f"  Output: {f}")
