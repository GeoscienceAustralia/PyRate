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
    params[cf.PARALLEL] = 0
    output_conf = Path(tmpdir).joinpath('conf.cfg')
    cf.write_config_file(params=params, output_conf_file=output_conf)
    params = configuration.Configuration(output_conf).__dict__

    dest_paths = [p.sampled_path for p in params[cf.INTERFEROGRAM_FILES]]
    # run conv2tif and prepifg, create the dest_paths files
    conv2tif.main(params)
    params[cf.INTERFEROGRAM_FILES].pop()
    prepifg.main(params)
    params[cf.INTERFEROGRAM_FILES].pop()
    preread_ifgs = process._create_ifg_dict(dest_paths, params=params)
    refpx, refpy = process._ref_pixel_calc(dest_paths, params)
    process._orb_fit_calc(params[cf.INTERFEROGRAM_FILES], params)
    process._ref_phase_estimation(dest_paths, params, refpx, refpy)

    maxvar, vcmt = process._maxvar_vcm_calc(dest_paths, params, preread_ifgs)

    # phase data after ref pixel has changed due to commit bf2f7ebd
    # Legacy tests won't match anymore
    np.testing.assert_array_almost_equal(maxvar, legacy_maxvar, decimal=4)
    np.testing.assert_array_almost_equal(legacy_vcm, vcmt, decimal=3)
    mpiops.run_once(shutil.rmtree, tmpdir)
