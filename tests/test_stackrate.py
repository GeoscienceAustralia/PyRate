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
# coding: utf-8
"""
This Python module contains tests for the stack.py PyRate module.
"""
import os
import shutil
import tempfile
import unittest

from numpy import eye, array, ones
import numpy as np
from numpy.testing import assert_array_almost_equal

import pyrate.core.orbital
import tests.common
from pyrate.core import shared, ref_phs_est as rpe, config as cf, covariance as vcm_module
from pyrate.core.stack import stack_rate
from pyrate import process, prepifg, conv2tif
from pyrate.configuration import Configuration
from tests.common import (SML_TEST_DIR, prepare_ifgs_without_phase,
    TEST_CONF_ROIPAC, pre_prepare_ifgs, remove_tifs)


def default_params():
    return {'pthr': 3, 'nsig': 3, 'maxsig': 2, 'parallel': 1, 'processes': 8}


class SinglePixelIfg(object):

    def __init__(self, timespan, phase):
        self.time_span = timespan
        self.phase_data = array([[phase]])


class StackRateTests(unittest.TestCase):
    """
    Tests the weighted least squares algorithm for determinining
    the best fitting velocity
    """

    def setUp(self):
        phase = [0.5, 3.5, 4, 2.5, 3.5, 1]
        timespan = [0.1, 0.7, 0.8, 0.5, 0.7, 0.2]
        self.ifgs = [SinglePixelIfg(s, p) for s, p in zip(timespan, phase)]

    def test_stack_rate(self):
        # Simple test with one pixel and equal weighting
        exprate = array([[5.0]])
        experr = array([[0.836242010007091]])
        expsamp = array([[5]])
        vcmt = eye(6, 6)
        mst = ones((6, 1, 1))
        mst[4] = 0
        params = default_params()
        rate, error, samples = stack_rate(self.ifgs, params, vcmt, mst)
        assert_array_almost_equal(rate, exprate)
        assert_array_almost_equal(error, experr)
        assert_array_almost_equal(samples, expsamp)


class LegacyEqualityTest(unittest.TestCase):
    """
    Tests equality with legacy data
    """

    @classmethod
    def setUpClass(cls):
        params = Configuration(TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = os.path.join(params[cf.OUT_DIR], cf.TMPDIR)
        shared.mkdir_p(params[cf.TMPDIR])
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, _, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST], params[cf.OBS_DIR])
        
        dest_paths = cf.get_dest_paths(base_ifg_paths, crop, params, xlks)
        # start run_pyrate copy
        ifgs = pre_prepare_ifgs(dest_paths, params)
        mst_grid = tests.common.mst_calculation(dest_paths, params)

        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)
        ifgs = prepare_ifgs_without_phase(dest_paths, params)
        for ifg in ifgs:
            ifg.close()
        _, ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)
        ifgs[0].open()
        r_dist = vcm_module.RDist(ifgs[0])()
        ifgs[0].close()
        maxvar = [vcm_module.cvd(i, params, r_dist)[0] for i in dest_paths]
        for ifg in ifgs:
            ifg.open()
        vcmt = vcm_module.get_vcmt(ifgs, maxvar)    
        for ifg in ifgs:
            ifg.close()     
            ifg.open()
        
        # Calculate stacked rate map
        params[cf.PARALLEL] = 1
        cls.rate, cls.error, cls.samples = tests.common.calculate_stack_rate(ifgs, params, vcmt, mst_mat=mst_grid)

        # Calculate stacked rate map
        params[cf.PARALLEL] = 0
        cls.rate_s, cls.error_s, cls.samples_s = tests.common.calculate_stack_rate(ifgs, params, vcmt, mst_mat=mst_grid)

        stackrate_dir = os.path.join(SML_TEST_DIR, 'stackrate')

        cls.rate_container = np.genfromtxt(os.path.join(stackrate_dir, 'stackmap.csv'), delimiter=',')
        cls.error_container = np.genfromtxt(os.path.join(stackrate_dir, 'errormap.csv'), delimiter=',')
        cls.samples_container = np.genfromtxt(os.path.join(stackrate_dir, 'coh_sta.csv'), delimiter=',')
    
        for ifg in ifgs:
            ifg.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)
        params = cf.get_config_params(TEST_CONF_ROIPAC)
        remove_tifs(params[cf.OBS_DIR])

    def test_stack_rate_full_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(self.rate, self.rate_s, decimal=3)

    def test_stackrate_error_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(self.error, self.error_s, decimal=3)

    def test_stackrate_samples_parallel(self):
        """
        python multiprocessing by rows vs serial
        """
        np.testing.assert_array_almost_equal(self.samples, self.samples_s, decimal=3)

    def test_stack_rate(self):
        """
        Compare with legacy data
        """
        np.testing.assert_array_almost_equal(self.rate_s, self.rate_container, decimal=3)

    def test_stackrate_error(self):
        """
        Compare with legacy data
        """
        np.testing.assert_array_almost_equal(self.error_s, self.error_container, decimal=3)

    def test_stackrate_samples(self):
        """
        Compare with legacy data
        """
        np.testing.assert_array_almost_equal(self.samples_s, self.samples_container, decimal=3)
