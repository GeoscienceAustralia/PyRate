#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
This Python module contains fixtures for use in the PyRate test suite.
"""
import shutil
import random
import string
import tempfile
import pytest

import pyrate.constants as C
from pyrate.constants import PYRATEPATH
from pyrate.core import mpiops, shared
from pyrate.configuration import Configuration
from tests.common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA, SML_TEST_DEM_TIF, MEXICO_CROPA_CONF
from tests.common import ROIPAC_SYSTEM_CONF, GAMMA_SYSTEM_CONF, GEOTIF_SYSTEM_CONF, SML_TEST_COH_LIST


@pytest.fixture
def tempdir():
    """
    tempdir for tests
    """
    def tmpdir():
        return tempfile.mkdtemp()
    return tmpdir


@pytest.fixture(params=[ROIPAC_SYSTEM_CONF, GEOTIF_SYSTEM_CONF, GAMMA_SYSTEM_CONF])
def system_conf(request):
    params = Configuration(request.param).__dict__
    yield request.param
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture
def random_filename(tmp_path_factory):
    def make_random_filename(ext=''):
        fname = ''.join(random.choice(string.ascii_lowercase) for _ in range(10))
        return tmp_path_factory.mktemp(basename='pyrate-').joinpath(fname + ext)
    return make_random_filename


# Make sure all MPI tests use this fixure
@pytest.fixture()
def mpisync(request):
    mpiops.comm.barrier()

    def fin():
        mpiops.comm.barrier()

    request.addfinalizer(fin)
    return mpiops.comm


@pytest.fixture(params=[0, 1])
def coh_mask(request):
    return request.param


@pytest.fixture(params=[1, 2])
def ref_est_method(request):
    return request.param


@pytest.fixture(params=[(-1, -1), (150.941666654, -34.218333314)])
def ref_pixel(request):
    return request.param


@pytest.fixture(params=[1, 3])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=C.ORB_METHOD_NAMES.keys())
def orbfit_method(request):
    return request.param


@pytest.fixture(params=C.ORB_DEGREE_NAMES.keys())
def orbfit_degrees(request):
    return request.param


@pytest.fixture(params=[1, 3])
def get_lks(request):
    return request.param


@pytest.fixture(params=[1, 4])
def get_crop(request):
    return request.param


@pytest.fixture()
def get_config():
    def params(conf_file):
        params_ = Configuration(conf_file).__dict__
        return params_
    return params


@pytest.fixture
def gamma_params():
    params = Configuration(TEST_CONF_GAMMA).__dict__
    shared.mkdir_p(params[C.OUT_DIR])
    yield params
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture
def roipac_params():
    params = Configuration(TEST_CONF_ROIPAC).__dict__
    yield params
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture
def mexico_cropa_params():
    params = Configuration(MEXICO_CROPA_CONF).__dict__
    yield params
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture(params=[TEST_CONF_GAMMA, TEST_CONF_ROIPAC])
def roipac_or_gamma_conf(request):
    return request.param


@pytest.fixture(params=[TEST_CONF_GAMMA], scope='session')
def gamma_conf(request):
    params = Configuration(TEST_CONF_GAMMA).__dict__
    yield request.param
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture
def coh_list_file():
    return SML_TEST_COH_LIST


@pytest.fixture
def dem():
    d = shared.dem_or_ifg(SML_TEST_DEM_TIF)
    d.open()
    return d


@pytest.fixture(params=[TEST_CONF_GAMMA, MEXICO_CROPA_CONF], scope='session')
def gamma_or_mexicoa_conf(request):
    params = Configuration(request.param).__dict__
    yield request.param
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)


@pytest.fixture(params=range(5))
def run_number(request):
    return request.param


GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs


@pytest.fixture
def ten_geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs[:10]


@pytest.fixture
def cropa_geotifs(mexico_cropa_params):
    m_paths = mexico_cropa_params[C.INTERFEROGRAM_FILES]
    return [m.converted_path for m in m_paths]
