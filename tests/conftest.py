import os
import random
import shutil
from pathlib import Path
import string
import tempfile
import pytest
from pyrate.core import config as cf, mpiops
from pyrate.configuration import Configuration
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
        params = cf.get_config_params(conf_file)
        # hack the rest of the params from Configuration class
        conf = Configuration(conf_file)
        params[cf.INTERFEROGRAM_FILES] = conf.interferogram_files
        params['rows'] = conf.rows
        params['cols'] = conf.cols
        return params
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


@pytest.fixture
def modify_config(tempdir, get_config, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees):
    def deliver_params(conf_file):
        params = get_config(conf_file)
        params[cf.IFG_LKSX] = get_lks
        params[cf.IFG_LKSY] = get_lks
        params[cf.IFG_CROP_OPT] = get_crop
        params[cf.ORBITAL_FIT_LOOKS_X] = orbfit_lks
        params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.IFG_FILE_LIST] = os.path.join(params[cf.OBS_DIR], 'interferogram_list.txt')
        tdir = tempdir()
        common.copytree(params[cf.OBS_DIR], tdir)
        params[cf.OBS_DIR] = tdir
        params[cf.OUT_DIR] = Path(tdir).joinpath('out')
        params[cf.PARALLEL] = True
        params[cf.APS_CORRECTION] = 0
        return params
    return deliver_params



@pytest.fixture(params=range(1, 2))
def get_lks(request):
    return request.param


@pytest.fixture(params=range(1, 2))
def get_crop(request):
    return request.param


@pytest.fixture(params=[1, 2])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=cf.ORB_METHOD_NAMES.values())
def orbfit_method(request):
    return request.param


@pytest.fixture(params=cf.ORB_DEGREE_NAMES.values())
def orbfit_degrees(request):
    return request.param
