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
import json
import click
from pyrate import pyratelog as pylog
from pyrate import config as cf
from pyrate.scripts import run_prepifg, run_pyrate, postprocessing
from pyrate import __version__

log = logging.getLogger(__name__)


@click.group()
@click.version_option(__version__, u'-V', u'--version')
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    """Commandline options and logging setup"""
    pylog.configure(verbosity)


@cli.command()
@click.argument('config_file')
def prepifg(config_file):
    """
    Convert input files to geotiff and perform multilooking
    (resampling) and/or cropping
    """
    config_file = os.path.abspath(config_file)
    params = cf.get_config_params(config_file)
    log.info('This job was run with the following parameters:')
    log.info(json.dumps(params, indent=4, sort_keys=True))
    run_prepifg.main(params)


@cli.command()
@click.argument('config_file')
@click.option('-r', '--rows', type=int, default=1,
              help='divide ifgs into this many rows')
@click.option('-c', '--cols', type=int, default=1,
              help='divide ifgs into this many columns')
def linrate(config_file, rows, cols):
    """
    Main PyRate workflow including time series and linear rate computation
    """
    config_file = os.path.abspath(config_file)
    _, dest_paths, params = cf.get_ifg_paths(config_file)
    log.info('This job was run with the following parameters:')
    log.info(json.dumps(params, indent=4, sort_keys=True))
    run_pyrate.process_ifgs(sorted(dest_paths), params, rows, cols)


@cli.command()
@click.argument('config_file')
@click.option('-r', '--rows', type=int, default=1,
              help='divide ifgs into this many rows. Must be same as '
                   'number of rows used previously in main workflow')
@click.option('-c', '--cols', type=int, default=1,
              help='divide ifgs into this many columns. Must be same as '
                   'number of cols used previously in main workflow')
def postprocess(config_file, rows, cols):
    """
    Reassemble PyRate output tiles and save as geotiffs
    """
    config_file = os.path.abspath(config_file)
    postprocessing.main(config_file, rows, cols)
