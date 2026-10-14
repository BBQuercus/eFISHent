"""Command line argument definitions and parsing."""

import argparse
import shutil
import sys
from typing import Any

import luigi

from . import __version__
from .constants import CLI_SHORTFORM, CONFIG_CLASSES
from .util import UniCode
from .validation import (
    existing_count_table,
    existing_fasta_file,
    existing_gtf_file,
    non_negative_int,
    percentage,
    positive_int,
    string_to_bool,
    validate_args,
)

GROUP_DESCRIPTIONS = {
    "General": "General configuration that will be used for all tasks.",
    "Run": "Options that change the behavior of the workflow.",
    "Sequence": "Details about the sequences the probe design will be performed on.",
    "Probe": "Probe filtering and design options.",
}
REQUIRED_PARAMS = ["reference_genome"]


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
