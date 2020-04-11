import os
import random
import shutil
import string
import tempfile
import pytest
from pyrate.core import config as cf, mpiops
from tests import common


@pytest.fixture
def tempdir():
    """
    tempdir for tests
    """
    def tmpdir():
        return tempfile.mkdtemp()
    return tmpdir


@pytest.fixture
def random_filename(tmpdir_factory):
    def make_random_filename(ext=''):
        dir = str(tmpdir_factory.mktemp('pyrate').realpath())
        fname = ''.join(random.choice(string.ascii_lowercase) for _ in range(10))
        return os.path.join(dir, fname + ext)
    return make_random_filename


@pytest.fixture()
def get_config():
    """
    Parameters
    ----------
    conf_file: str
        config file

    Returns
    -------
    params: dict
        dict of params
    """
    def params(conf_file):
        return cf.get_config_params(conf_file)
    return params


# Make sure all MPI tests use this fixure
@pytest.fixture()
def mpisync(request):
    mpiops.comm.barrier()

    def fin():
        mpiops.comm.barrier()

    request.addfinalizer(fin)
    return mpiops.comm


@pytest.fixture(params=[0, 1])
def roipac_or_gamma(request):
    return request.param


@pytest.fixture(params=[1, 2])
def ref_est_method(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def row_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def col_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def modify_config(request, tempdir, get_config):
    test_conf = common.TEST_CONF_ROIPAC
    params_dict = get_config(test_conf)
    params_dict[cf.IFG_LKSX] = request.param
    params_dict[cf.IFG_LKSY] = request.param
    params_dict[cf.OBS_DIR] = tempdir()
    common.copytree(common.SML_TEST_GAMMA, params_dict[cf.OBS_DIR])
    params_dict[cf.IFG_FILE_LIST] = os.path.join(params_dict[cf.OBS_DIR], 'ifms_17')
    params_dict[cf.PARALLEL] = False
    params_dict[cf.APS_CORRECTION] = 0
    yield params_dict
    # clean up
    shutil.rmtree(params_dict[cf.OBS_DIR])


@pytest.fixture(params=range(1, 6))
def get_lks(request):
    return request.param


@pytest.fixture(params=range(1, 3))
def get_crop(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=[cf.INDEPENDENT_METHOD, cf.NETWORK_METHOD])
def orbfit_method(request):
    return request.param


@pytest.fixture(params=[cf.PLANAR, cf.QUADRATIC, cf.PART_CUBIC])
def orbfit_degrees(request):
    return request.param
