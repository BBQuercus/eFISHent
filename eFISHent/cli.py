import argparse
import configparser
import logging
import os
import sys
import tempfile

import luigi

from .alignment import BuildBowtieIndex
from .cleanup import CleanUpOutput
from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig
from .kmers import BuildJellyfishIndex
from .util import UniCode

CONFIG_CLASSES = [GeneralConfig, RunConfig, SequenceConfig, ProbeConfig]

GROUP_DESCRIPTIONS = {
    f"{UniCode.blue} General": "General configuration that will be used for all tasks.",
    f"{UniCode.green} Run": "Options that change the behavior of the workflow.",
    f"{UniCode.red} Sequence": "Details about the sequences the probe design will be performed on.",
    f"{UniCode.magenta} Probe": "Probe filtering and design options.",
}


def string_to_bool(value):
    """Workaround for using typed boolean values as arguments."""
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif value.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def get_parameter_type(param: luigi.Parameter) -> type:
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
    required_params = ["reference_genome"]

    for group, config_class in zip(groups, CONFIG_CLASSES):
        for name, param in config_class().get_params():
            group.add_argument(
                f"--{name.replace('_', '-')}",
                type=get_parameter_type(param),
                required=name in required_params,
                default=param._default,
                help=f"{param.description} "
                f"[default: {'-' if not param._default else param._default}, "
                f"required: {name in required_params}]",
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
    config.read("luigi.cfg")

    for section, config_class in zip(
        ["GeneralConfig", "RunConfig", "SequenceConfig", "ProbeConfig"], CONFIG_CLASSES
    ):
        for name in config_class().get_param_names():
            config.set(section, name, str(vars(args).get(name)))
    with open(config_file, "w") as f:
        config.write(f)


def set_logging_level(verbose: bool, debug: bool) -> logging.Logger:
    """Set the logging level of luigi and custom logger."""
    LOG_FORMAT = "%(asctime)s %(levelname)-4s [%(name)s] %(message)s"
    logging.basicConfig(format=LOG_FORMAT, force=True)

    if debug:
        luigi_level = "DEBUG"
        custom_level = logging.DEBUG
    elif verbose:
        luigi_level = "WARNING"
        custom_level = logging.INFO
    else:
        luigi_level = "WARNING"
        custom_level = logging.WARNING

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

        tasks = []
        if args.build_indices:
            tasks = [BuildJellyfishIndex(), BuildBowtieIndex()]
        else:
            tasks = [CleanUpOutput()]

        luigi.build(tasks, local_scheduler=True)

    if tasks[-1].complete():
        logger.info(f"{UniCode.party} eFISHent has finished running!")
