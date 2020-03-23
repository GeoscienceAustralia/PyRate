# coding: utf-8
#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
import glob
import os
import random
import shutil
import string
import tempfile
from os.path import join, basename

import numpy as np
import pytest

from . import common
import conv2tif
import core.orbital
import core.shared
import prepifg
import process
from common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA
from core import mpiops, config as cf
from configuration import Configuration

TRAVIS = True if "TRAVIS" in os.environ else False

legacy_maxvar = [
    15.4156637191772,
    2.85829424858093,
    34.3486289978027,
    2.59190344810486,
    3.18510007858276,
    3.61054635047913,
    1.64398515224457,
    14.9226036071777,
    5.13451862335205,
    6.82901763916016,
    10.9644861221313,
    14.5026779174805,
    29.3710079193115,
    8.00364685058594,
    2.06328082084656,
    5.66661834716797,
    5.62802362442017,
]


@pytest.fixture()
def tempdir():
    """tempdir for tests"""

    def tmpdir():
        """ """
        return tempfile.mkdtemp()

    return tmpdir


@pytest.fixture
def random_filename(tmpdir_factory):
    """

    Args:
      tmpdir_factory: 

    Returns:

    """

    def make_random_filename(ext=""):
        """

        Args:
          ext: Default value = "")

        Returns:

        """
        dir = str(tmpdir_factory.mktemp("pyrate").realpath())
        fname = "".join(random.choice(string.ascii_lowercase) for _ in range(10))
        return os.path.join(dir, fname + ext)

    return make_random_filename


@pytest.fixture()
def get_config():
    """

    Args:

    Returns:
      dict: params** -- dict of params

    """

    def params(conf_file):
        """

        Args:
          conf_file: 

        Returns:

        """

        return Configuration(conf_file).__dict__

    return params


# Make sure all MPI tests use this fixure
@pytest.fixture()
def mpisync(request):
    """

    Args:
      request: 

    Returns:

    """
    mpiops.comm.barrier()

    def fin():
        """ """
        mpiops.comm.barrier()

    request.addfinalizer(fin)
    return mpiops.comm


@pytest.fixture(params=[0, 1])
def roipac_or_gamma(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[1, 2])
def ref_est_method(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[1, 2, 5])
def row_splits(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[1, 2, 5])
def col_splits(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[1, 2, 5])
def modify_config(request, tempdir, get_config):
    """

    Args:
      request: param tempdir:
      get_config: 
      tempdir: 

    Returns:

    """
    test_conf = common.TEST_CONF_ROIPAC
    params_dict = get_config(test_conf)
    params_dict[cf.IFG_LKSX] = request.param
    params_dict[cf.IFG_LKSY] = request.param
    params_dict[cf.OUT_DIR] = tempdir()
    common.copytree(common.SML_TEST_GAMMA, params_dict[cf.OUT_DIR])
    params_dict[cf.IFG_FILE_LIST] = os.path.join(params_dict[cf.OUT_DIR], "ifms_17")
    params_dict[cf.PARALLEL] = False
    params_dict[cf.APS_CORRECTION] = 0
    yield params_dict
    # clean up
    shutil.rmtree(params_dict[cf.OUT_DIR])


@pytest.fixture(params=range(1, 6))
def get_lks(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=range(1, 3))
def get_crop(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


def test_vcm_legacy_vs_mpi(mpisync, tempdir, get_config):
    from common import SML_TEST_DIR, TEST_CONF_ROIPAC

    params_dict = get_config(TEST_CONF_ROIPAC)
    LEGACY_VCM_DIR = os.path.join(SML_TEST_DIR, 'vcm')
    legacy_vcm = np.genfromtxt(os.path.join(LEGACY_VCM_DIR, 'vcmt.csv'), delimiter=',')
    if mpiops.rank == 0:
        outdir = tempdir()
    else:
        outdir = None
    outdir = mpiops.comm.bcast(outdir, root=0)
    params_dict[cf.OUT_DIR] = outdir
    params_dict[cf.PARALLEL] = False
    xlks, ylks, crop = cf.transform_params(params_dict)

    dest_paths = []
    for interferogram_file in params_dict["interferogram_files"]:
        dest_paths.append(interferogram_file.sampled_path)

    # run prepifg, create the dest_paths files
    if mpiops.rank == 0:
        conv2tif.main(params_dict)
        prepifg.main(params_dict)

    mpiops.comm.barrier()

    tiles = core.shared.get_tiles(dest_paths[0], rows=1, cols=1)
    preread_ifgs = process._create_ifg_dict(dest_paths, params=params_dict, tiles=tiles)
    refpx, refpy = process._ref_pixel_calc(dest_paths, params_dict)
    process._orb_fit_calc(dest_paths, params_dict)
    process._ref_phase_estimation(dest_paths, params_dict, refpx, refpy)

    maxvar, vcmt = process._maxvar_vcm_calc(dest_paths, params_dict, preread_ifgs)
    np.testing.assert_array_almost_equal(maxvar, legacy_maxvar, decimal=4)
    # np.testing.assert_array_almost_equal(legacy_vcm, vcmt, decimal=3)
    if mpiops.rank == 0:
        shutil.rmtree(params_dict[cf.OUT_DIR])


@pytest.fixture(params=[1, 2, 5])
def orbfit_lks(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[cf.INDEPENDENT_METHOD, cf.NETWORK_METHOD])
def orbfit_method(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


@pytest.fixture(params=[cf.PLANAR, cf.QUADRATIC, cf.PART_CUBIC])
def orbfit_degrees(request):
    """

    Args:
      request: 

    Returns:

    """
    return request.param


def _tifs_same(dir1, dir2, tif):
    """

    Args:
      dir1: param dir2:
      tif: 
      dir2: 

    Returns:

    """
    stack_tif_s = os.path.join(dir1, tif)
    stack_tif_m = os.path.join(dir2, tif)
    common.assert_ifg_phase_equal(stack_tif_m, stack_tif_s)


def test_prepifg_mpi(mpisync, get_config, tempdir, roipac_or_gamma, get_lks, get_crop):
    """

    Args:
      mpisync: param get_config:
      tempdir: param roipac_or_gamma:
      get_lks: param get_crop:
      get_config: 
      roipac_or_gamma: 
      get_crop: 

    Returns:

    """

    if roipac_or_gamma == 1:
        params = get_config(TEST_CONF_GAMMA)
    else:
        params = get_config(TEST_CONF_ROIPAC)
    outdir = mpiops.run_once(tempdir)
    params[cf.OUT_DIR] = outdir
    params[cf.PARALLEL] = False
    params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
    params[cf.IFG_CROP_OPT] = get_crop
    conv2tif.main(params)
    prepifg.main(params)
    common.remove_tifs(params[cf.OUT_DIR])

    if mpiops.rank == 0:
        if roipac_or_gamma == 1:
            params_s = get_config(TEST_CONF_GAMMA)
        else:
            params_s = get_config(TEST_CONF_ROIPAC)
        params_s[cf.OUT_DIR] = tempdir()
        params_s[cf.PARALLEL] = True
        params_s[cf.IFG_LKSX], params_s[cf.IFG_LKSY] = get_lks, get_lks
        params_s[cf.IFG_CROP_OPT] = get_crop
        conv2tif.main(params)
        if roipac_or_gamma == 1:
            base_unw_paths = glob.glob(join(common.SML_TEST_GAMMA, "*_utm.unw"))
            prepifg.main(params)
        else:
            base_unw_paths = glob.glob(join(common.SML_TEST_OBS, "*.unw"))
            prepifg.main(params_s)

        mpi_tifs = glob.glob(join(outdir, "*.tif"))
        serial_tifs = glob.glob(join(params[cf.OUT_DIR], "*.tif"))
        mpi_tifs.sort()
        serial_tifs.sort()
        # 17 geotifs, and 17 mlooked tifs
        assert len(mpi_tifs) == len(serial_tifs)
        for m_f, s_f in zip(mpi_tifs, serial_tifs):
            assert basename(m_f) == basename(s_f)

        shutil.rmtree(outdir)
        shutil.rmtree(params_s[cf.OUT_DIR])
