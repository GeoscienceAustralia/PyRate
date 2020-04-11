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
from pyrate.main import conv2tif_handler, prepifg_handler, process_handler, merge_handler
from pyrate.constants import PYRATEPATH
from pyrate.core.config import write_config_file, get_config_params
from tests.common import ROIPAC_SYSTEM_CONF, GAMMA_SYSTEM_CONF, GEOTIF_SYSTEM_CONF


# @pytest.fixture
# def roipac_conf():
#     input_config_path =
#     params =

def test_roipac_workflow():
    # joblib_conf_path =  input_config_path.
    conv2tif_handler(ROIPAC_SYSTEM_CONF)
    prepifg_handler(ROIPAC_SYSTEM_CONF)
    process_handler(ROIPAC_SYSTEM_CONF)
    merge_handler(ROIPAC_SYSTEM_CONF)


def test_gamma_workflow():

    conv2tif_handler(GAMMA_SYSTEM_CONF)
    prepifg_handler(GAMMA_SYSTEM_CONF)
    process_handler(GAMMA_SYSTEM_CONF)
    merge_handler(GAMMA_SYSTEM_CONF)


def test_geotiff_workflow():
    prepifg_handler(GEOTIF_SYSTEM_CONF)
    process_handler(GEOTIF_SYSTEM_CONF)
    merge_handler(GEOTIF_SYSTEM_CONF)
