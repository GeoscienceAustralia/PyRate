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
import glob
import shutil
import numpy as np
import os
import pytest

import pyrate.core.orbital
import pyrate.core.shared
from pyrate import process, prepifg, conv2tif
from pyrate.core import mpiops, config as cf
from pyrate.core.shared import Ifg
from tests import common
from tests.common import SML_TEST_DIR, TEST_CONF_ROIPAC, TEST_CONF_GAMMA
from tests.test_covariance import legacy_maxvar


@pytest.fixture()
def get_mpi_config():
    def params(conf_file):
        params = cf.get_config_params(conf_file)
        return params
    return params


@pytest.fixture(params=[TEST_CONF_ROIPAC, TEST_CONF_GAMMA])
def roipac_or_gamma(request):
    return request.param


@pytest.fixture(params=range(1, 6))
def get_lks(request):
    return request.param


@pytest.fixture(params=range(1, 3))
def get_crop(request):
    return request.param


def test_vcm_legacy_vs_mpi(mpisync, tempdir, get_mpi_config):
    params_dict = get_mpi_config(TEST_CONF_ROIPAC)
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

    # run conv2tif and prepifg, create the dest_paths files
    conv2tif.main(params_dict)
    prepifg.main(params_dict)

    tiles = pyrate.core.shared.get_tiles(dest_paths[0], rows=1, cols=1)
    preread_ifgs = process._create_ifg_dict(dest_paths, params=params_dict, tiles=tiles)
    refpx, refpy = process._ref_pixel_calc(dest_paths, params_dict)
    process._orb_fit_calc(dest_paths, params_dict)
    process._ref_phase_estimation(dest_paths, params_dict, refpx, refpy)

    maxvar, vcmt = process._maxvar_vcm_calc(dest_paths, params_dict, preread_ifgs)
    np.testing.assert_array_almost_equal(maxvar, legacy_maxvar, decimal=4)
    np.testing.assert_array_almost_equal(legacy_vcm, vcmt, decimal=3)
    mpiops.run_once(shutil.rmtree, outdir)
    mpiops.run_once(common.remove_tifs, params_dict[cf.OBS_DIR])


def _tifs_same(dir1, dir2, tif):
    stackrate_tif_s = os.path.join(dir1, tif)
    stackrate_tif_m = os.path.join(dir2, tif)
    common.assert_ifg_phase_equal(stackrate_tif_m, stackrate_tif_s)


def test_conv2tif_prepifg_mpi(mpisync, get_mpi_config, tempdir, roipac_or_gamma, get_lks, get_crop):
    from os.path import join, basename
    params = get_mpi_config(TEST_CONF_GAMMA)
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
    mpiops.comm.barrier()
    mpiops.run_once(common.remove_tifs, params[cf.OBS_DIR])

    if mpiops.rank == 0:
        params_s = get_mpi_config(TEST_CONF_GAMMA)
        params_s[cf.OUT_DIR] = tempdir()
        params_s[cf.PARALLEL] = True
        params_s[cf.IFG_LKSX], params_s[cf.IFG_LKSY] = get_lks, get_lks
        params_s[cf.IFG_CROP_OPT] = get_crop
        conv2tif.main(params)
        if roipac_or_gamma == 1:
            prepifg.main(params)
        else:
            prepifg.main(params_s)

        mpi_tifs = glob.glob(join(outdir, "*.tif"))
        serial_tifs = glob.glob(join(params[cf.OUT_DIR], "*.tif"))
        mpi_tifs.sort()
        serial_tifs.sort()
        # 17 geotifs, and 17 mlooked tifs + 1 dem
        assert len(mpi_tifs) == len(serial_tifs) == 18
        for m_f, s_f in zip(mpi_tifs, serial_tifs):
            assert basename(m_f) == basename(s_f)
            common.assert_tifs_equal(m_f, s_f)

        shutil.rmtree(outdir)
        shutil.rmtree(params_s[cf.OUT_DIR])
        common.remove_tifs(params[cf.OBS_DIR])


# TODO: mpi test for the rest of the pipeline
def test_convert2tif_mpi():
    pass


def test_process_mpi():
    pass


def test_merge_mpi():
    pass
