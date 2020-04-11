import os
import random
import string
import tempfile
import pytest
from pyrate.core import mpiops


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


@pytest.fixture(params=[1, 2, 5])
def row_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def col_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def orbfit_lks(request):
    return request.param
