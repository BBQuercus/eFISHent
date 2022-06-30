"""Command line interface."""

from pathlib import Path
import argparse
import configparser
import logging
import os
import sys
import tempfile
from typing import Any, List

import luigi

from .alignment import BuildBowtieIndex
from .cleanup import CleanUpOutput
from .constants import CONFIG_CLASSES
from .kmers import BuildJellyfishIndex
from .util import UniCode


GROUP_DESCRIPTIONS = {
    f"{UniCode.blue} General": "General configuration that will be used for all tasks.",
    f"{UniCode.green} Run": "Options that change the behavior of the workflow.",
    f"{UniCode.red} Sequence": "Details about the sequences the probe design will be performed on.",
    f"{UniCode.magenta} Probe": "Probe filtering and design options.",
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

    raise argparse.ArgumentTypeError("Boolean value expected.")


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
    utility = parser.add_argument_group(f"{UniCode.cyan}General utilities{UniCode.end}")
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
        version="%(prog)s 0.0.1",
        help="Show %(prog)s's version number.",
    )
    utility.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Set program output to verbose printing all important steps.",
    )
    utility.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)


def _add_groups(parser: argparse.ArgumentParser) -> None:
    """Add the main option groups to the parser."""
    groups = [
        parser.add_argument_group(
            f"{name} options{UniCode.end}", description=description
        )
        for name, description in GROUP_DESCRIPTIONS.items()
    ]

    for group, config_class in zip(groups, CONFIG_CLASSES):
        for name, param in config_class().get_params():
            param_type = get_parameter_type(param)
            is_required = name in REQUIRED_PARAMS
            default = (
                "-"
                if (not param._default and param_type != string_to_bool)
                else param._default
            )
            group.add_argument(
                f"--{name.replace('_', '-')}",
                type=param_type,
                required=is_required,
                default=param._default,
                help=f"{param.description} "
                f"[default: {default}, "
                f"required: {is_required}]",
            )


def _parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog="eFISHent",
        description=f"{UniCode.bold}eFISHent {UniCode.fishing} {UniCode.dna} to design all your probes.{UniCode.end}",
        epilog=(
            'See the online wiki at "https://github.com/BBQuercus/eFISHent/wiki" for an overview.\n'
            f"We hope you enjoy using eFISHent {UniCode.party}!"
        ),
        add_help=False,
    )
    _add_groups(parser)
    _add_utilities(parser)
    try:
        if len(sys.argv) == 1:
            parser.print_help()
            parser.exit(0)
    except Exception as e:
        print(e)
    return parser.parse_args()


def create_custom_config(args: argparse.Namespace, config_file: str) -> None:
    """Create a custom config file."""
    config = configparser.ConfigParser()
    config_path = Path(__file__).resolve().parent.joinpath("luigi.cfg").as_posix()
    config.read(config_path)

    for section, config_class in zip(
        ["GeneralConfig", "RunConfig", "SequenceConfig", "ProbeConfig"], CONFIG_CLASSES
    ):
        for name in config_class().get_param_names():
            value = vars(args).get(name)
            config.set(section, name, str(value))
            if name == "threads":
                threads = min(vars(args).get(name), os.cpu_count())  # type: ignore
                config.set(section, name, str(threads))

    with open(config_file, "w") as f:
        config.write(f)


def set_logging_level(verbose: bool, debug: bool) -> logging.Logger:
    """Set the logging level of luigi and custom logger."""
    log_format = "%(asctime)s %(levelname)-4s %(message)s"

    if debug:
        log_format = (
            "%(asctime)s %(levelname)-4s [%(name)s] "
            "%(filename)s %(funcName)s %(lineno)d / %(thread)d - %(message)s"
        )
        luigi_level = "DEBUG"
        custom_level = logging.DEBUG
        logfile = "eFISHent.log"
    elif verbose:
        luigi_level = "WARNING"
        custom_level = logging.INFO
        logfile = None
    else:
        luigi_level = "WARNING"
        custom_level = logging.WARNING
        logfile = None

    logging.basicConfig(filename=logfile, format=log_format, force=True)
    logging.getLogger("luigi").setLevel(luigi_level)
    logging.getLogger("luigi-interface").setLevel(luigi_level)
    luigi.interface.core.log_level = luigi_level

    logger = logging.getLogger("custom-logger")
    logger.setLevel(custom_level)
    return logger


def main():
    """Run eFISHent tasks."""
    args = _parse_args()
    logger = set_logging_level(args.verbose, args.debug)
    logger.info(f"{UniCode.fishing} eFISHent has started running...")

    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        create_custom_config(args, config_file)
        luigi.configuration.add_config_path(config_file)

        tasks: List[luigi.Task] = []
        if args.build_indices:
            tasks = [BuildJellyfishIndex(), BuildBowtieIndex()]
        else:
            tasks = [CleanUpOutput()]

        luigi.build(tasks, local_scheduler=True)

    if tasks[-1].complete():
        logger.info(f"{UniCode.party} eFISHent has finished running!")
