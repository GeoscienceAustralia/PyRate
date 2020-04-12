# coding: utf-8
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
This Python module contains tests for mpi operations in PyRate.
Tun this module as 'mpirun -n 4 pytest tests/test_mpi.py'
"""
import pytest
import shutil
from pathlib import Path

from pyrate.main import conv2tif_handler, prepifg_handler, process_handler, merge_handler
from pyrate.core.config import OBS_DIR, PROCESSOR
from pyrate.core import config as cf
from pyrate import conv2tif, prepifg, process, merge
from pyrate.configuration import Configuration

from tests.common import ROIPAC_SYSTEM_CONF, GAMMA_SYSTEM_CONF, GEOTIF_SYSTEM_CONF, copytree


@pytest.fixture(params=[2])
def get_lks(request):
    return request.param


@pytest.fixture(params=[2])
def get_crop(request):
    return request.param


@pytest.fixture(params=[2])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=cf.ORB_METHOD_NAMES.keys())
def orbfit_method(request):
    return request.param


@pytest.fixture(params=cf.ORB_DEGREE_NAMES.keys())
def orbfit_degrees(request):
    return request.param


@pytest.fixture(params=[ROIPAC_SYSTEM_CONF, GEOTIF_SYSTEM_CONF, GAMMA_SYSTEM_CONF])
def conf(request):
    return request.param


@pytest.fixture
def get_params(get_config, conf):
    return get_config(conf)


@pytest.fixture()
def get_config(tempdir, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees):
    def modify_params(conf_file):
        params = cf.get_config_params(conf_file)

        tdir = Path(tempdir())
        copytree(params[cf.OBS_DIR], tdir)

        # manipulate params
        params[cf.OBS_DIR] = tdir.as_posix()
        params[cf.OUT_DIR] = tdir.joinpath('out')
        params[cf.PARALLEL] = 1
        params[cf.APS_CORRECTION] = 0
        params[cf.IFG_LKSX] = get_lks
        params[cf.IFG_LKSY] = get_lks
        params[cf.IFG_CROP_OPT] = get_crop
        params[cf.ORBITAL_FIT_LOOKS_X] = orbfit_lks
        params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees

        # write new temp config
        output_conf = tdir.joinpath('config.conf')
        cf.write_config_file(params=params, output_conf_file=output_conf)

        # hack the rest of the params from Configuration class
        conf = Configuration(output_conf)
        return conf.__dict__
    return modify_params


def test_workflow(get_params):
    params = get_params
    if params[PROCESSOR] != 2:  # if not geotif
        conv2tif.main(params)
    prepifg.main(params)
    process.main(params)
    merge.main(params)
    shutil.rmtree(params[OBS_DIR])

# TODO:  MPI comparison system tests


def test_gamma_workflow():
    """check the handlers are working as expected"""
    conv2tif_handler(GAMMA_SYSTEM_CONF)
    prepifg_handler(GAMMA_SYSTEM_CONF)
    process_handler(GAMMA_SYSTEM_CONF)
    merge_handler(GAMMA_SYSTEM_CONF)
