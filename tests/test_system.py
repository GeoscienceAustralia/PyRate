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
from pyrate.main import conv2tif_handler, prepifg_handler, process_handler, merge_handler
from pyrate.core.config import OBS_DIR, PROCESSOR
from pyrate import conv2tif, prepifg, process, merge

from tests.common import ROIPAC_SYSTEM_CONF, GAMMA_SYSTEM_CONF, GEOTIF_SYSTEM_CONF


@pytest.fixture(params=[GEOTIF_SYSTEM_CONF, GAMMA_SYSTEM_CONF, ROIPAC_SYSTEM_CONF])
def conf(request):
    return request.param


@pytest.fixture
def roipac_or_gamma_params(modify_config, conf):
    return modify_config(conf)


def test_workflow(roipac_or_gamma_params):
    params_dict = roipac_or_gamma_params
    if params_dict[PROCESSOR] != 2:  # if not geotif
        conv2tif.main(params_dict)
    prepifg.main(params_dict)
    process.main(params_dict)
    merge.main(params_dict)
    shutil.rmtree(params_dict[OBS_DIR])


def test_gamma_workflow():
    """One test to wrap the handlers"""
    conv2tif_handler(GAMMA_SYSTEM_CONF)
    prepifg_handler(GAMMA_SYSTEM_CONF)
    process_handler(GAMMA_SYSTEM_CONF)
    merge_handler(GAMMA_SYSTEM_CONF)
