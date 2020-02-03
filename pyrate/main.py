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
import logging
import argparse
from argparse import RawTextHelpFormatter

from constants import CLI_DESCRIPTION, CONV2TIF, PREPIFG, PROCESS, MERGE
from core import config as cf
import conv2tif, prepifg, process, merge
from core import pyratelog
from core import user_experience
import time
import multiprocessing
from shutil import copyfile

from core.user_experience import break_number_into_factors
from core.config import OBS_DIR, OUT_DIR
from configuration import Configuration
import pathlib

# Turn off MPI warning
os.environ['OMPI_MCA_btl_base_warn_component_unused'] = '0'
import logging

log = logging.getLogger(__name__)

log.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
# fh = logging.FileHandler('logfile.log')
# fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
# log.addHandler(fh)
log.addHandler(ch)


def conv2tif_handler(config_file):
    """
    Convert interferograms to geotiff.
    """
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    conv2tif.main(config.__dict__)


def prepifg_handler(config_file):
    """
    Perform multilooking and cropping on geotiffs.
    """
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    prepifg.main(config.__dict__)


def process_handler(config_file):
    """
    Time series and linear rate computation.
    """
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    process.main(config.__dict__)


def merge_handler(config_file):
    """
    Reassemble computed tiles and save as geotiffs.
    """
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    merge.main(config.__dict__)
    user_experience.delete_tsincr_files(config.__dict__)


def main():

    start_time = time.time()
    log.debug("Starting PyRate")

    parser = argparse.ArgumentParser(prog='pyrate', description=CLI_DESCRIPTION, add_help=True, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--verbosity', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help="Increase output verbosity")

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    parser_conv2tif = subparsers.add_parser('conv2tif', help='Convert interferograms to geotiff.', add_help=True)
    parser_conv2tif.add_argument('-f', '--config_file', action="store", type=str, default=None,  help="Pass configuration file", required=True)

    parser_prepifg = subparsers.add_parser('prepifg', help='Perform multilooking and cropping on geotiffs.', add_help=True)
    parser_prepifg.add_argument('-f', '--config_file', action="store", type=str, default=None, help="Pass configuration file", required=True)

    parser_process = subparsers.add_parser('process', help='Main processing workflow including corrections, time series and stacking computation.', add_help=True)
    parser_process.add_argument('-f', '--config_file', action="store", type=str, default=None, help="Pass configuration file", required=True)

    parser_merge = subparsers.add_parser('merge', help="Reassemble computed tiles and save as geotiffs.", add_help=True)
    parser_merge.add_argument('-f', '--config_file', action="store", type=str, default=None, help="Pass configuration file", required=False)

    parser_workflow = subparsers.add_parser('workflow', help="Run all the PyRate processes", add_help=True)
    parser_workflow.add_argument('-f', '--config_file', action="store", type=str, default=None, help="Pass configuration file", required=False)

    args = parser.parse_args()

    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    if args.verbosity:
        log.setLevel(args.verbosity)
        log.info("Verbosity set to " + str(args.verbosity) + ".")

    if args.command == "conv2tif":
        conv2tif_handler(args.config_file)

    if args.command == "prepifg":
        prepifg_handler(args.config_file)

    if args.command == "process":
        process_handler(args.config_file)

    if args.command == "merge":
        merge_handler(args.config_file)

    if args.command == "workflow":
        log.info("***********CONV2TIF**************")
        conv2tif_handler(args.config_file)

        log.info("***********PREPIFG**************")
        prepifg_handler(args.config_file)

        log.info("***********PROCESS**************")
        process_handler(args.config_file)

        log.info("***********MERGE**************")
        merge_handler(args.config_file)

    log.debug("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
