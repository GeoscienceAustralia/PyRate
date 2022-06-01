#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module defines executable run configuration for the PyRate software
"""

import os
import argparse
from argparse import RawTextHelpFormatter
import time
from pathlib import Path

import pyrate.constants as C
from pyrate.constants import CLI_DESCRIPTION
from pyrate import conv2tif, prepifg, correct, merge
from pyrate.core.logger import pyratelogger as log, configure_stage_log
from pyrate.core import mpiops
from pyrate.configuration import Configuration
from pyrate.core.shared import mpi_vs_multiprocess_logging
from pyrate.core.stack import stack_calc_wrapper
from pyrate.core.timeseries import timeseries_calc_wrapper


def _params_from_conf(config_file):
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    params = config.__dict__
    return params


def update_params_due_to_ifg_selection(config):
    params = config.__dict__
    if config.phase_closure_filtered_ifgs_list(params).exists():
        params = config.refresh_ifg_list(params)
        correct._create_ifg_dict(params)
        correct._update_params_with_tiles(params)
    return params


# pylint: disable=too-many-statements
# JUSTIFICATION: extracting statements into separate functions would make it harder to understand.
def main():
    """The main PyRate runner program, refer to `docs/index.rst` and `--help`"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='pyrate', description=CLI_DESCRIPTION, add_help=True,
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        '-v', '--verbosity', type=str,
        default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help="Increase output verbosity"
    )

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    parser_conv2tif = subparsers.add_parser(
        'conv2tif',
        help='<Optional> Convert interferograms to geotiff.',
        add_help=True
    )

    parser_prepifg = subparsers.add_parser(
        'prepifg',
        help='Perform multilooking, cropping and coherence masking to interferogram geotiffs.',
        add_help=True
    )

    parser_correct = subparsers.add_parser(
        'correct',
        help='Calculate and apply corrections to interferogram phase data.',
        add_help=True
    )

    parser_ts = subparsers.add_parser(
        'timeseries',
        help='<Optional> Timeseries inversion of interferogram phase data.',
        add_help=True
    )

    parser_stack = subparsers.add_parser(
        'stack',
        help='<Optional> Stacking of interferogram phase data.',
        add_help=True
    )

    parser_merge = subparsers.add_parser(
        'merge',
        help="Reassemble computed tiles and save as geotiffs.",
        add_help=True
    )

    parser_workflow = subparsers.add_parser(
        'workflow',
        help="<Optional> Sequentially run all the PyRate processing steps.",
        add_help=True
    )

    sub_parsers = [
        parser_conv2tif,
        parser_prepifg,
        parser_correct,
        parser_merge,
        parser_ts,
        parser_stack,
        parser_workflow
    ]

    for sub_parser in sub_parsers:
        sub_parser.add_argument('-f', '--config_file', action="store", type=str, default=None,
                       help="Pass configuration file", required=False)

    args = parser.parse_args()

    params = mpiops.run_once(_params_from_conf, args.config_file)

    log_path = Path(params[C.OUT_DIR]).joinpath('pyrate.log.').as_posix()
    configure_stage_log(args.verbosity, args.command, log_path)

    log.debug("Starting PyRate")
    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    if args.verbosity:
        log.setLevel(args.verbosity)
        log.info(f"Verbosity set to {args.verbosity}.")

    if args.command == "conv2tif":
        conv2tif.main(params)

    if args.command == "prepifg":
        prepifg.main(params)

    if args.command == "correct":
        config_file = os.path.abspath(args.config_file)
        config = Configuration(config_file)
        correct.main(config)

    if args.command == "timeseries":
        config_file = os.path.abspath(args.config_file)
        config = Configuration(config_file)
        timeseries(config)

    if args.command == "stack":
        config_file = os.path.abspath(args.config_file)
        config = Configuration(config_file)
        stack(config)

    if args.command == "merge":
        merge.main(params)

    if args.command == "workflow":
        log.info("***********CONV2TIF**************")
        conv2tif.main(params)

        log.info("***********PREPIFG**************")
        params = mpiops.run_once(_params_from_conf, args.config_file)
        prepifg.main(params)

        log.info("***********CORRECT**************")
        # reset params as prepifg modifies params
        config_file = os.path.abspath(args.config_file)
        config = Configuration(config_file)
        correct.main(config)

        log.info("***********TIMESERIES**************")
        config = Configuration(config_file)
        timeseries(config)

        log.info("***********STACK**************")
        config = Configuration(config_file)
        stack(config)

        log.info("***********MERGE**************")
        params = mpiops.run_once(_params_from_conf, args.config_file)
        merge.main(params)

    elapsed_secs = time.time() - start_time
    log.info(f"--- Runtime = {elapsed_secs} seconds ---")


def timeseries(config: Configuration) -> None:
    """The runner command for calculating the timeseries of a file set"""
    params = config.__dict__
    mpi_vs_multiprocess_logging("timeseries", params)
    params = update_params_due_to_ifg_selection(config=config)
    timeseries_calc_wrapper(params)


def stack(config: Configuration) -> None:
    """The runner command for stacking a set of files"""
    params = config.__dict__
    mpi_vs_multiprocess_logging("stack", params)
    params = update_params_due_to_ifg_selection(config=config)
    stack_calc_wrapper(params)


if __name__ == "__main__":
    main()
