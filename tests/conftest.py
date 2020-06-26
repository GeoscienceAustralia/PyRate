import shutil
import random
import string
import tempfile
import pytest
from pyrate.core import mpiops, config as cf, prepifg_helper, shared
from pyrate.configuration import Configuration
from tests.common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA
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
    shutil.rmtree(params[cf.OUT_DIR])


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


@pytest.fixture(params=cf.ORB_METHOD_NAMES.keys())
def orbfit_method(request):
    return request.param


@pytest.fixture(params=cf.ORB_DEGREE_NAMES.keys())
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
        params = cf.get_config_params(conf_file)
        return params
    return params


@pytest.fixture
def gamma_params():
    params = Configuration(TEST_CONF_GAMMA).__dict__
    shutil.rmtree(params[cf.OUT_DIR])
    shared.mkdir_p(params[cf.OUT_DIR])
    yield params
    shutil.rmtree(params[cf.OUT_DIR])


@pytest.fixture
def roipac_params():
    params = Configuration(TEST_CONF_ROIPAC).__dict__
    yield params
    shutil.rmtree(params[cf.OUT_DIR])


@pytest.fixture(params=[TEST_CONF_GAMMA, TEST_CONF_ROIPAC])
def roipac_or_gamma_conf(request):
    return request.param


@pytest.fixture(params=[TEST_CONF_GAMMA], scope='session')
def gamma_conf(request):
    params = Configuration(TEST_CONF_GAMMA).__dict__
    yield request.param
    shutil.rmtree(params[cf.OUT_DIR], ignore_errors=True)


@pytest.fixture
def coh_list_file():
    return SML_TEST_COH_LIST
