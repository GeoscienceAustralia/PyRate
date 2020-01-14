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
import shutil
import numpy as np
import pytest
import os
import tempfile
import random
import string

import core.orbital
import core.shared
import common
import process, prepifg, merge, conv2tif
from common import (small_data_setup, reconstruct_mst,  reconstruct_stack_rate, SML_TEST_DEM_HDR_GAMMA, pre_prepare_ifgs)
import common
from test_covariance import legacy_maxvar
from core import algorithm, ref_phs_est as rpe, mpiops, config as cf, covariance, refpixel

TRAVIS = True if 'TRAVIS' in os.environ else False


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
    base_unw_paths = cf.original_ifg_paths(params_dict[cf.IFG_FILE_LIST], params_dict[cf.OBS_DIR])
    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params_dict, xlks)

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
    np.testing.assert_array_almost_equal(legacy_vcm, vcmt, decimal=3)
    if mpiops.rank == 0:
        shutil.rmtree(outdir)
        common.remove_tifs(params_dict[cf.OBS_DIR])


@pytest.fixture(params=[1, 2, 5])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=[cf.INDEPENDENT_METHOD, cf.NETWORK_METHOD])
def orbfit_method(request):
    return request.param


@pytest.fixture(params=[cf.PLANAR, cf.QUADRATIC, cf.PART_CUBIC])
def orbfit_degrees(request):
    return request.param


def _tifs_same(dir1, dir2, tif):
    stack_tif_s = os.path.join(dir1, tif)
    stack_tif_m = os.path.join(dir2, tif)
    common.assert_ifg_phase_equal(stack_tif_m, stack_tif_s)

def test_prepifg_mpi(mpisync, get_config, tempdir, roipac_or_gamma, get_lks, get_crop):
    from common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA
    from os.path import join, basename
    if roipac_or_gamma == 1:
        params = get_config(TEST_CONF_GAMMA)
    else:
        params = get_config(TEST_CONF_ROIPAC)
    outdir = mpiops.run_once(tempdir)
    params[cf.OUT_DIR] = outdir
    params[cf.PARALLEL] = False
    params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
    params[cf.IFG_CROP_OPT] = get_crop
    if roipac_or_gamma == 1:
        params[cf.IFG_FILE_LIST] = join(common.SML_TEST_GAMMA, 'ifms_17')
        params[cf.OBS_DIR] = common.SML_TEST_GAMMA
        params[cf.DEM_FILE] = common.SML_TEST_DEM_GAMMA
        params[cf.DEM_HEADER_FILE] = common.SML_TEST_DEM_HDR_GAMMA
    conv2tif.main(params)
    prepifg.main(params)
    common.remove_tifs(params[cf.OBS_DIR])    

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
        common.remove_tifs(params[cf.OBS_DIR])
