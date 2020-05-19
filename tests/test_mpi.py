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
import shutil
import numpy as np
import os
from pathlib import Path

import pyrate.core.orbital
import pyrate.core.shared
from pyrate import process, prepifg, conv2tif, configuration
from pyrate.core import mpiops, config as cf
from tests import common
from tests.common import SML_TEST_DIR
from tests.test_covariance import legacy_maxvar


def test_vcm_legacy_vs_mpi(mpisync, tempdir, roipac_or_gamma_conf):
    params = configuration.Configuration(roipac_or_gamma_conf).__dict__
    LEGACY_VCM_DIR = os.path.join(SML_TEST_DIR, 'vcm')
    legacy_vcm = np.genfromtxt(os.path.join(LEGACY_VCM_DIR, 'vcmt.csv'), delimiter=',')
    tmpdir = Path(mpiops.run_once(tempdir))
    mpiops.run_once(common.copytree, params[cf.OBS_DIR], tmpdir)
    params[cf.OUT_DIR] = tmpdir.joinpath('out')
    params[cf.PARALLEL] = False
    xlks, ylks, crop = cf.transform_params(params)
    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST], params[cf.OBS_DIR])
    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    # run conv2tif and prepifg, create the dest_paths files
    conv2tif.main(params)
    prepifg.main(params)

    tiles = pyrate.core.shared.get_tiles(dest_paths[0], rows=1, cols=1)
    preread_ifgs = process._create_ifg_dict(dest_paths, params=params, tiles=tiles)
    refpx, refpy = process._ref_pixel_calc(dest_paths, params)
    process._orb_fit_calc(dest_paths, params)
    process._ref_phase_estimation(dest_paths, params, refpx, refpy)

    maxvar, vcmt = process._maxvar_vcm_calc(dest_paths, params, preread_ifgs)
    np.testing.assert_array_almost_equal(maxvar, legacy_maxvar, decimal=4)
    np.testing.assert_array_almost_equal(legacy_vcm, vcmt, decimal=3)
    mpiops.run_once(shutil.rmtree, tmpdir)
