"""Command line interface."""

from pathlib import Path
from typing import Any, List, TYPE_CHECKING
import argparse
import configparser
import logging
import os
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
    from .indexing import BuildBowtieIndex
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
    "sequence_file": existing_fasta_file,
    "encode_count_table": existing_count_table,
    "analyze_probeset": existing_fasta_file,
    # Positive integers
    "threads": positive_int,
    "min_length": positive_int,
    "max_length": positive_int,
    "kmer_length": positive_int,
    "max_kmers": positive_int,
    "optimization_time_limit": positive_int,
    # Non-negative integers
    "spacing": non_negative_int,
    "max_off_targets": non_negative_int,
    "sequence_similarity": non_negative_int,
    # Percentages
    "min_gc": percentage,
    "max_gc": percentage,
    "formamide_concentration": percentage,
    "max_expression_percentage": percentage,
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
    # Boolean parameters - show yes/no
    "build_indices": "yes/no",
    "save_intermediates": "yes/no",
    "is_plus_strand": "yes/no",
    "is_endogenous": "yes/no",
    "no_alternative_loci": "yes/no",
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
    utility.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)


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

    if errors:
        parser.error("\n  " + "\n  ".join(errors))


def _parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    import shutil

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

    args = parser.parse_args()
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
    from .console import pipeline_progress, print_completion, print_header

    start_time = time.time()
    args = _parse_args()
    logger = set_logging_level(args.silent, args.debug)

    # Print header banner
    if not args.silent:
        print_header(__version__)
    else:
        logger.info(f"{UniCode.fishing} eFISHent v{__version__} starting...")

    # Lazy imports - only load heavy modules after arg parsing
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtieIndex
    from .kmers import BuildJellyfishIndex

    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        create_custom_config(args, config_file)
        luigi.configuration.add_config_path(config_file)

        tasks: List[luigi.Task] = []
        if args.build_indices:
            tasks = [BuildJellyfishIndex(), BuildBowtieIndex()]
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
                with pipeline_progress(total_stages=8):
                    luigi.build(tasks, local_scheduler=True)
        else:
            luigi.build(tasks, local_scheduler=True)

    duration = time.time() - start_time

    if tasks[-1].complete():
        if not args.silent:
            print_completion(format_duration(duration))
        else:
            logger.info(
                f"{UniCode.party} eFISHent has finished running in {format_duration(duration)}!"
            )
