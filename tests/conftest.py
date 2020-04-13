import os
import random
import string
import tempfile
import pytest
from pyrate.core import mpiops, config as cf
from tests.common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA
from tests.common import ROIPAC_SYSTEM_CONF, GAMMA_SYSTEM_CONF, GEOTIF_SYSTEM_CONF


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
    return request.param


@pytest.fixture
def random_filename(tmpdir_factory):
    def make_random_filename(ext=''):
        dir = str(tmpdir_factory.mktemp('pyrate').realpath())
        fname = ''.join(random.choice(string.ascii_lowercase) for _ in range(10))
        return os.path.join(dir, fname + ext)
    return make_random_filename


# Make sure all MPI tests use this fixure
@pytest.fixture()
def mpisync(request):
    mpiops.comm.barrier()

    def fin():
        mpiops.comm.barrier()

    request.addfinalizer(fin)
    return mpiops.comm


@pytest.fixture(params=[1, 2])
def ref_est_method(request):
    return request.param


@pytest.fixture(params=[1, 3])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=cf.ORB_METHOD_NAMES.keys())
def orbfit_method(request):
    return request.param


@pytest.fixture(params=cf.ORB_DEGREE_NAMES.keys())
def orbfit_degrees(request):
    return request.param


@pytest.fixture(params=[1, 3])
def get_lks(request):
    return request.param


@pytest.fixture(params=[1, 3])
def get_crop(request):
    return request.param


@pytest.fixture()
def get_config():
    def params(conf_file):
        params = cf.get_config_params(conf_file)
        return params
    return params


@pytest.fixture(params=[TEST_CONF_GAMMA, TEST_CONF_ROIPAC])
def roipac_or_gamma_conf(request):
    return request.param
