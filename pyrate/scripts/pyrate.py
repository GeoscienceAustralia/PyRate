import click
import logging
from os.path import abspath
from pyrate import pyratelog as pylog
from pyrate import config as cf
from pyrate.scripts import run_prepifg, run_pyrate

log = logging.getLogger(__name__)


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    pylog.configure(verbosity)


@cli.command()
@click.argument('config_file')
def prepifg(config_file):
    config_file = abspath(config_file)
    params = cf.get_config_params(config_file)
    run_prepifg.main(params)


@cli.command()
@click.argument('config_file')
def linrate(config_file):
    config_file = abspath(config_file)
    # params = cf.get_config_params(config_file)
    run_pyrate.main(config_file)
