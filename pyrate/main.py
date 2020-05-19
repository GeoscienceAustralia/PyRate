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

from pyrate.constants import CLI_DESCRIPTION
from pyrate import conv2tif, prepifg, process, merge
from pyrate.core.logger import pyratelogger as log, configure_stage_log
from pyrate.core import config as cf
from pyrate.core import mpiops
from pyrate.configuration import Configuration


def _params_from_conf(config_file):
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    return config.__dict__


def main():

    start_time = time.time()

    parser = argparse.ArgumentParser(prog='pyrate', description=CLI_DESCRIPTION, add_help=True,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--verbosity', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help="Increase output verbosity")

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    parser_conv2tif = subparsers.add_parser('conv2tif', help='Convert interferograms to geotiff.', add_help=True)
    parser_conv2tif.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                 help="Pass configuration file", required=True)

    parser_prepifg = subparsers.add_parser('prepifg', help='Perform multilooking and cropping on geotiffs.',
                                           add_help=True)
    parser_prepifg.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                help="Pass configuration file", required=True)

    parser_process = subparsers.add_parser(
        'process', help='Main processing workflow including corrections, time series and stacking computation.',
        add_help=True)
    parser_process.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                help="Pass configuration file", required=True)

    parser_merge = subparsers.add_parser('merge', help="Reassemble computed tiles and save as geotiffs.",
                                         add_help=True)
    parser_merge.add_argument('-f', '--config_file', action="store", type=str, default=None,
                              help="Pass configuration file", required=False)

    parser_workflow = subparsers.add_parser('workflow', help="Run all the PyRate processes", add_help=True)
    parser_workflow.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                 help="Pass configuration file", required=False)

    args = parser.parse_args()

    params = mpiops.run_once(_params_from_conf, args.config_file)

    configure_stage_log(args.verbosity, args.command, Path(params[cf.OUT_DIR]).joinpath('pyrate.log.').as_posix())

    log.debug("Starting PyRate")
    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    if args.verbosity:
        log.setLevel(args.verbosity)
        log.info("Verbosity set to " + str(args.verbosity) + ".")

    if args.command == "conv2tif":
        conv2tif.main(params)

    if args.command == "prepifg":
        prepifg.main(params)

    if args.command == "process":
        process.main(params)

    if args.command == "merge":
        merge.main(params)

    if args.command == "workflow":
        log.info("***********CONV2TIF**************")
        conv2tif.main(params)

        log.info("***********PREPIFG**************")
        params = mpiops.run_once(_params_from_conf, args.config_file)
        prepifg.main(params)

        log.info("***********PROCESS**************")
        # reset params as prepifg modifies params
        params = mpiops.run_once(_params_from_conf, args.config_file)
        process.main(params)

        # process might modify params too
        params = mpiops.run_once(_params_from_conf, args.config_file)
        log.info("***********MERGE**************")
        merge.main(params)

    log.debug("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
