import pytest
import os
import tempfile
import random
import string
from pyrate import config as cf
from pyrate import mpiops


@pytest.fixture()
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
        fname = ''.join(random.choice(string.ascii_lowercase)
                        for _ in range(10))
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


@pytest.fixture(params=[1, 5, 10])
def row_splits(request):
    return request.param


@pytest.fixture(params=[1, 5, 10])
def col_splits(request):
    return request.param
