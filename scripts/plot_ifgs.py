#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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


def main():

    start_time = time.time()

    parser = argparse.ArgumentParser(prog='pyrate', description=CLI_DESCRIPTION, add_help=True,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--verbosity', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help="Increase output verbosity")

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    parser_plot = subparsers.add_parser('plot_ifgs', help='Plot interferogram.', add_help=True)
    parser_plot.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                help="Pass configuration file", required=True)

    args = parser.parse_args()

    params = mpiops.run_once(_params_from_conf, args.config_file)

    configure_stage_log(args.verbosity, args.command, Path(params[C.OUT_DIR]).joinpath('pyrate.log.').as_posix())

    log.debug("Plotting")
    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    log.info("--- Runtime = %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
