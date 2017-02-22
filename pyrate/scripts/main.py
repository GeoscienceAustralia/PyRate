"""Pyrate execution main module"""

from os.path import abspath
import logging
import click
from pyrate import pyratelog as pylog
from pyrate import config as cf
from pyrate.scripts import run_prepifg, run_pyrate, postprocessing

log = logging.getLogger(__name__)


@click.group()
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
    Prepifg wrapper. This performs both geotif coverserion and
    also multilooking (resampling) and cropping
    """
    config_file = abspath(config_file)
    params = cf.get_config_params(config_file)
    run_prepifg.main(params)


@cli.command()
@click.argument('config_file')
@click.option('-r', '--rows', type=int, default=1,
              help='divide ifgs into this many rows')
@click.option('-c', '--cols', type=int, default=1,
              help='divide ifgs into this many columns')
def linrate(config_file, rows, cols):
    """ Time series and linear rate computation"""
    config_file = abspath(config_file)
    run_pyrate.main(config_file, rows, cols)


@cli.command()
@click.argument('config_file')
@click.option('-r', '--rows', type=int, default=1,
              help='divide ifgs into this many rows. Must be same as '
                   'number of rows used in linrate')
@click.option('-c', '--cols', type=int, default=1,
              help='divide ifgs into this many columns. Must be same as '
                   'number of cols used in linrate')
def postprocess(config_file, rows, cols):
    """ assemble the linrate and timeseries tiles back together"""
    config_file = abspath(config_file)
    postprocessing.main(config_file, rows, cols)
