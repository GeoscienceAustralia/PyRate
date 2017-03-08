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
This Python module contains tests for the linrate.py PyRate module.
"""
from __future__ import print_function
import os
import shutil
import sys
import tempfile
import unittest
from numpy import eye, array, ones

import numpy as np
from numpy.testing import assert_array_almost_equal

import pyrate.orbital
import tests.common
from pyrate import config as cf
from pyrate import ref_phs_est as rpe
from pyrate import shared
from pyrate import vcm as vcm_module
from pyrate.linrate import linear_rate
from pyrate.scripts import run_pyrate, run_prepifg
from tests.common import SYD_TEST_DIR, prepare_ifgs_without_phase


def default_params():
    return {'pthr': 3, 'nsig': 3, 'maxsig': 2, 'parallel': 1, 'processes': 8}


class SinglePixelIfg(object):

    def __init__(self, timespan, phase):
        self.time_span = timespan
        self.phase_data = array([[phase]])


class LinearRateTests(unittest.TestCase):
    """
    Tests the weighted least squares algorithm for determinining
    the best fitting velocity
    """

    def setUp(self):
        phase = [0.5, 3.5, 4, 2.5, 3.5, 1]
        timespan = [0.1, 0.7, 0.8, 0.5, 0.7, 0.2]
        self.ifgs = [SinglePixelIfg(s, p) for s, p in zip(timespan, phase)]

    def test_linear_rate(self):
        # Simple test with one pixel and equal weighting
        exprate = array([[5.0]])
        experr = array([[0.836242010007091]])  # from Matlab Pirate
        expsamp = array([[5]])
        vcmt = eye(6, 6)
        mst = ones((6, 1, 1))
        mst[4] = 0
        params = default_params()
        rate, error, samples = linear_rate(self.ifgs, params, vcmt, mst)
        assert_array_almost_equal(rate, exprate)
        assert_array_almost_equal(error, experr)
        assert_array_almost_equal(samples, expsamp)


class MatlabEqualityTest(unittest.TestCase):
    """
    Tests equality vs matlab
    """

    @classmethod
    def setUpClass(cls):
        params = cf.get_config_params(os.path.join(SYD_TEST_DIR,
                                                   'pyrate_system_test.conf'))

        cls.temp_out_dir = tempfile.mkdtemp()

        sys.argv = ['run_prepifg.py',
                    os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf')]
        params[cf.OUT_DIR] = cls.temp_out_dir
        run_prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, _, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(
            params[cf.IFG_FILE_LIST])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop, params, xlks)

        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = tests.common.mst_calculation(dest_paths, params)

        refx, refy = run_pyrate.ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.orbital.remove_orbital_error(ifgs, params)
        ifgs = prepare_ifgs_without_phase(dest_paths, params)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        maxvar = [vcm_module.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm_module.get_vcmt(ifgs, maxvar)

        # Calculate linear rate map
        params[cf.PARALLEL] = 1
        cls.rate, cls.error, cls.samples = tests.common.calculate_linear_rate(
            ifgs, params, vcmt, mst_mat=mst_grid)

        params[cf.PARALLEL] = 2
        cls.rate_2, cls.error_2, cls.samples_2 = \
            tests.common.calculate_linear_rate(ifgs, params, vcmt,
                                               mst_mat=mst_grid)

        params[cf.PARALLEL] = 0
        # Calculate linear rate map
        cls.rate_s, cls.error_s, cls.samples_s = \
            tests.common.calculate_linear_rate(ifgs, params, vcmt,
                                               mst_mat=mst_grid)

        matlab_linrate_dir = os.path.join(SYD_TEST_DIR, 'matlab_linrate')

        cls.rate_matlab = np.genfromtxt(
            os.path.join(matlab_linrate_dir, 'stackmap.csv'), delimiter=',')
        cls.error_matlab = np.genfromtxt(
            os.path.join(matlab_linrate_dir, 'errormap.csv'), delimiter=',')

        cls.samples_matlab = np.genfromtxt(
            os.path.join(matlab_linrate_dir, 'coh_sta.csv'), delimiter=',')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_linear_rate_full_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(
            self.rate, self.rate_s, decimal=3)

    def test_linrate_error_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(
            self.error, self.error_s, decimal=3)

    def test_linrate_samples_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(
            self.samples, self.samples_s, decimal=3)

    def test_linrate_full_parallel_pixel(self):
        """
        python multiprocessing by pixel vs serial
        """
        np.testing.assert_array_almost_equal(
            self.rate_2, self.rate_s, decimal=3)

    def test_linrate_error_parallel_pixel(self):
        """
        python multiprocessing by pixel vs serial
        """
        np.testing.assert_array_almost_equal(
            self.error_2, self.error_s, decimal=3)

    def test_linrate_samples_parallel_pixel(self):
        """
        python multiprocessing pixel level vs serial
        """
        np.testing.assert_array_almost_equal(
            self.samples_2, self.samples_s, decimal=3)

    def test_linear_rate(self):
        """
        python vs matlab
        """
        np.testing.assert_array_almost_equal(
            self.rate_s, self.rate_matlab, decimal=3)

    def test_linrate_error(self):
        """
        python vs matlab
        """
        np.testing.assert_array_almost_equal(
            self.error_s, self.error_matlab, decimal=3)

    def test_linrate_samples(self):
        """
        python lin rate samples vs matlab
        """
        np.testing.assert_array_almost_equal(
            self.samples_s, self.samples_matlab, decimal=3)


if __name__ == "__main__":
    unittest.main()
