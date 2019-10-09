#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
from pyrate.core import config as cf
from pyrate import (converttogeotif, prepifg, process, postprocess)
from pyrate import CONV2TIF, PREPIFG, PROCESS, MERGE # Step names
from pyrate.core import pyratelog
import time

log = logging.getLogger(__name__)

def converttogeotiff_handler(config_file):
    """
    Convert interferograms to geotiff.
    """
    config_file = os.path.abspath(config_file)
    params = cf.get_config_params(config_file, step=CONV2TIF)
    converttogeotif.main(params)


def prepifg_handler(config_file):
    """
    Perform multilooking and cropping on geotiffs.
    """
    config_file = os.path.abspath(config_file)
    params = cf.get_config_params(config_file, step=PREPIFG)
    prepifg.main(params)


def process_handler(config_file, rows, cols):
    """
    Time series and linear rate computation.
    """
    config_file = os.path.abspath(config_file)
    _, dest_paths, params = cf.get_ifg_paths(config_file, step=PROCESS)
    process.process_ifgs(sorted(dest_paths), params, rows, cols)


def postprocess_handler(config_file, rows, cols):
    """
    Reassemble computed tiles and save as geotiffs.
    """
    config_file = os.path.abspath(config_file)
    _, _, params = cf.get_ifg_paths(config_file, step=MERGE)
    postprocess.main(params, rows, cols)

CLI_DESC = """
PyRate workflow: 

    Step 1: converttogeotiff
    Step 2: prepifg
    Step 3: process
    Step 4: postprocess 

Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for 
more details.
"""

def main():
    start_time = time.time()
    log.debug("Starting PyRate")
    parser = argparse.ArgumentParser(prog='pyrate',
                                     description=CLI_DESC,
                                     add_help=True, 
                                     formatter_class=RawTextHelpFormatter)

    # Have to add verbosity argument to each subparser so define args here
    #  to make life easier.
    verbosity_args = ['-v', '--verbosity']
    verbosity_kwargs = {'type': str,
                        'default': 'INFO',
                        'choices': ['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        'help': "Increase output verbosity"}
    
    parser.add_argument(*verbosity_args, **verbosity_kwargs)

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    # create the parser for the "converttogeotiff" command
    parser_converttogeotiff = \
        subparsers.add_parser('conv2tif',
                              help='Convert interferograms to geotiff.', 
                              add_help=True)

    parser_converttogeotiff.add_argument('-f', '--config_file', action="store", 
                                         type=str, default=None, 
                                         help="Pass configuration file", 
                                         required=True)

    parser_converttogeotiff.add_argument(*verbosity_args, **verbosity_kwargs)

    # create the parser for the "prepifg" command
    parser_prepifg = subparsers.add_parser('prepifg',
                                           help=('Perform multilooking and '
                                                 'cropping on geotiffs.'), 
                                           add_help=True)

    parser_prepifg.add_argument('-f', '--config_file', action="store", 
                                type=str, default=None, 
                                help="Pass configuration file", required=True)

    parser_prepifg.add_argument(*verbosity_args, **verbosity_kwargs)

    # create the parser for the "process" command
    parser_process = subparsers.add_parser('process',
                                           help=('Main processing workflow '
                                                 'including corrections, time '
                                                 'series and linear rate '
                                                 'computation.'),
                                           add_help=True)

    parser_process.add_argument('-f', '--config_file', action="store", 
                                type=str, default=None,
                                help="Pass configuration file", 
                                required=True)

    parser_process.add_argument('-r', '--rows', type=int, required=False, default=1,
                                help=("divide ifgs into this many rows. Must "
                                      "be same as number of rows used "
                                      "previously in main workflow."))

    parser_process.add_argument('-c', '--cols', type=int, required=False, default=1,
                                help=("divide ifgs into this many columns. "
                                      "Must be same as number of cols used "
                                      "previously in main workflow."))
        
    parser_process.add_argument(*verbosity_args, **verbosity_kwargs)
 
    # create the parser for the "postprocess" command
    parser_postprocess = subparsers.add_parser('merge',
                                           help=("Reassemble computed tiles "
                                                 "and save as geotiffs."), 
                                           add_help=True)

    parser_postprocess.add_argument('-f', '--config_file', action="store", 
                                type=str, default=None,
                                help="Pass configuration file", required=False)

    parser_postprocess.add_argument('-r', '--rows', type=int, required=False, default=1,
                                help=("divide ifgs into this many rows. Must "
                                      "be same as number of rows used "
                                      "previously in main workflow."))

    parser_postprocess.add_argument('-c', '--cols', type=int, required=False, default=1,
                                help=("divide ifgs into this many columns. "
                                      "Must be same as number of cols used "
                                      "previously in main workflow."))

    parser_postprocess.add_argument(*verbosity_args, **verbosity_kwargs)

    args = parser.parse_args()

    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    if args.verbosity:
        pyratelog.configure(args.verbosity)
        log.info("Verbosity set to " + str(args.verbosity) + ".")

    if args.command == "conv2tif":
        converttogeotiff_handler(args.config_file)

    if args.command == "prepifg":
        prepifg_handler(args.config_file)

    if args.command == "process":
        process_handler(args.config_file, args.rows, args.cols)

    if args.command == "merge":
        postprocess_handler(args.config_file, args.rows, args.cols)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    main()

