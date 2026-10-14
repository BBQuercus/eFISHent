"""Command line interface."""

import argparse
import configparser
import logging
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

# Silence some troublemakers before heavy imports
warnings.filterwarnings("ignore")
logging.getLogger("luigi-interface").setLevel(level=logging.CRITICAL)
logging.getLogger("matplotlib.font_manager").setLevel(level=logging.CRITICAL)

import luigi

from . import __version__
from .arguments import _parse_args
from .constants import CONFIG_CLASSES
from .dependencies import (
    _ensure_deps_on_path,
    check_all_dependencies,
    check_required_dependencies,
)
from .parameter_warnings import validate_parameter_warnings
from .runner import _run_tasks, _task_classes
from .updates import check_for_updates, self_update
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

if TYPE_CHECKING:
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
    from .indexing import BuildBowtie2Index, BuildBowtieIndex
    from .kmers import BuildJellyfishIndex


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
