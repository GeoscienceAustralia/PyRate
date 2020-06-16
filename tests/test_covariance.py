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
This Python module contains tests for the covariance.py PyRate module.
"""
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from numpy import array
import numpy as np
from numpy.testing import assert_array_almost_equal

from pyrate.core import shared, ref_phs_est as rpe, ifgconstants as ifc, config as cf
from pyrate import process, prepifg, conv2tif
from pyrate.core.covariance import cvd, get_vcmt, RDist
from pyrate.configuration import Configuration
import pyrate.core.orbital
from pyrate.core import roipac
from pyrate.core.config import parse_namelist
from tests import common
from tests.common import (small5_mock_ifgs, small5_ifgs, TEST_CONF_ROIPAC,
    small_data_setup, prepare_ifgs_without_phase)


class CovarianceTests(unittest.TestCase):
    def setUp(self):
        self.ifgs = small_data_setup()
        for i in self.ifgs:
            i.mm_converted = True
        params = dict()
        params[cf.NO_DATA_VALUE] = 0
        params[cf.NAN_CONVERSION] = True
        self.params = params
        self.r_dist = RDist(self.ifgs[0])()

    def test_covariance_basic(self):
        ifgs = small5_ifgs()
        for i in ifgs:
            i.open()

            if bool((i.phase_data == 0).all()) is True:
                raise Exception("All zero")

            maxvar, alpha = cvd(i, self.params, self.r_dist, calc_alpha=True)
            self.assertTrue(maxvar is not None)
            self.assertTrue(alpha is not None)
            print("maxvar: %s, alpha: %s" % (maxvar, alpha))

    def test_covariance_17ifgs(self):
        # After raw data import
        # (no reference pixel correction and units in radians)
        exp_maxvar = [5.6149, 8.7710, 2.9373, 0.3114, 12.9931, 2.0459, 0.4236,
                      2.1243, 0.4745, 0.6725, 0.8333, 3.8232, 3.3052, 2.4925,
                      16.0159, 2.8025, 1.4345]

        exp_alpha = [0.0356, 0.1601, 0.5128, 0.5736, 0.0691, 0.1337, 0.2333,
                    0.3202, 1.2338, 0.4273, 0.9024, 0.1280, 0.3585, 0.1599,
                    0.0110, 0.1287, 0.0676]

        act_maxvar = []
        act_alpha = []
        for i in self.ifgs:

            if bool((i.phase_data == 0).all()) is True:
                raise Exception("All zero")

            maxvar, alpha = cvd(i, self.params, self.r_dist, calc_alpha=True)
            self.assertTrue(maxvar is not None)
            self.assertTrue(alpha is not None)
           
            act_maxvar.append(maxvar)
            act_alpha.append(alpha)

        assert_array_almost_equal(act_maxvar, exp_maxvar, decimal=3)

        # This test fails for greater than 1 decimal place.
        # Discrepancies observed in distance calculations.
        assert_array_almost_equal(act_alpha, exp_alpha, decimal=1)


class VCMTests(unittest.TestCase):

    def setUp(self):
        self.ifgs = small_data_setup()

    def test_vcm_basic(self):
        ifgs = small5_mock_ifgs(5, 9)
        maxvar = [8.486, 12.925, 6.313, 0.788, 0.649]

        exp = array([[8.486, 5.2364, 0.0, 0.0, 0.0],
            [5.2364, 12.925,  4.5165,  1.5957,  0.0],
            [0.0, 4.5165, 6.313, 1.1152, 0.0],
            [0.0, 1.5957, 1.1152, 0.788, -0.3576],
            [0.0, 0.0, 0.0, -0.3576, 0.649]])

        act = get_vcmt(ifgs, maxvar)
        assert_array_almost_equal(act, exp, decimal=3)

    def test_vcm_17ifgs(self):
        # TODO: maxvar should be calculated by vcm.cvd
        maxvar = [2.879, 4.729, 22.891, 4.604, 3.290, 6.923, 2.519, 13.177,
                  7.548, 6.190, 12.565, 9.822, 18.484, 7.776, 2.734, 6.411,
                  4.754]

        exp = array([
            [2.879, 0.0, -4.059, -1.820, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 4.729, 0.0, 0.0, 1.972, 0.0, 0.0, -3.947, -2.987, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-4.059, 0.0, 22.891, 5.133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             -7.497, -10.285, 0.0, 0.0, 0.0, 0.0],
            [-1.820, 0.0, 5.133, 4.604, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             3.362, 0.0, 0.0, -1.774, 0.0, 0.0],
            [0.0, 1.972, 0.0, 0.0, 3.290, 2.386, 1.439, -3.292, -2.492, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 2.386, 6.923, 2.088, 0.0, 0.0, -3.273,
             -4.663, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.439, 2.088, 2.519, 0.0, 0.0, 1.974, 0.0,
             0.0, 0.0, -2.213, 0.0, 0.0, 0.0],
            [0.0, -3.947, 0.0, 0.0, -3.292, 0.0, 0.0, 13.177, 4.986, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 4.596, -3.957],
            [0.0, -2.987, 0.0, 0.0, -2.492, 0.0, 0.0, 4.986, 7.548, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 2.995],
            [0.0, 0.0, 0.0, 0.0, 0.0, -3.273, 1.974, 0.0, 0.0, 6.190, 4.410,
             0.0, 0.0, -3.469, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, -4.663, 0.0, 0.0, 0.0, 4.410, 12.565,
             0.0, 0.0, 4.942, 0.0, 0.0, 0.0],
            [0.0, 0.0, -7.497, 3.362, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             9.8221, 6.737, 0.0, -2.591, 0.0, 0.0],
            [0.0, 0.0, -10.285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.737,
             18.484, 0.0, 3.554, -5.443, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.213, 0.0, 0.0, -3.469, 4.942,
             0.0, 0.0, 7.776, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, -1.774, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.591,
             3.554, 0.0, 2.734, -2.093, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.596, 0.0, 0.0, 0.0, 0.0,
             -5.443, 0.0, -2.093, 6.411, -2.760],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.957, 2.995, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, -2.760, 4.754]
        ])

        act = get_vcmt(self.ifgs, maxvar)
        assert_array_almost_equal(act, exp, decimal=3)


legacy_maxvar = [15.4156637191772,
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
                 5.62802362442017]


class LegacyEqualityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        params = Configuration(TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = os.path.join(cls.temp_out_dir, cf.TMPDIR)
        shared.mkdir_p(params[cf.TMPDIR])
        params[cf.REF_EST_METHOD] = 2
        conv2tif.main(params)
        prepifg.main(params)
        cls.params = params
        base_ifg_paths = [c.unwrapped_path for c in params[cf.INTERFEROGRAM_FILES]]
        dest_paths = [c.converted_path for c in params[cf.INTERFEROGRAM_FILES]]
        dest_paths = dest_paths[:-2]
        for i in dest_paths:
            Path(i).chmod(0o664)  # assign write permission as conv2tif output is readonly
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        refx, refy = process._ref_pixel_calc(dest_paths, params)
        headers = [roipac.roipac_header(i, cls.params) for i in base_ifg_paths]
        pyrate.core.orbital.remove_orbital_error(ifgs, params, headers)
        ifgs = prepare_ifgs_without_phase(dest_paths, params)
        for ifg in ifgs:
            ifg.close()
        _, cls.ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)
        ifgs[0].open()
        r_dist = RDist(ifgs[0])()
        ifgs[0].close()
        # Calculate interferogram noise
        cls.maxvar = [cvd(i, params, r_dist, calc_alpha=True, save_acg=True, write_vals=True)[0] for i in dest_paths]
        cls.vcmt = get_vcmt(ifgs, cls.maxvar)
        for ifg in ifgs:
            ifg.close()

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        shutil.rmtree(cls.temp_out_dir)
        params = cf.get_config_params(TEST_CONF_ROIPAC)
        common.remove_tifs(params[cf.OBS_DIR])

    def test_legacy_maxvar_equality_small_test_files(self):
        np.testing.assert_array_almost_equal(self.maxvar, legacy_maxvar,
                                             decimal=3)

    def test_legacy_vcmt_equality_small_test_files(self):
        from tests.common import SML_TEST_DIR
        LEGACY_VCM_DIR = os.path.join(SML_TEST_DIR, 'vcm')
        legacy_vcm = np.genfromtxt(os.path.join(LEGACY_VCM_DIR, 'vcmt.csv'), delimiter=',')
        np.testing.assert_array_almost_equal(legacy_vcm, self.vcmt, decimal=3)

    def test_metadata(self):
        for ifg in self.ifgs:
            if not ifg.is_open:
                ifg.open()
            assert ifc.PYRATE_MAXVAR in ifg.meta_data
            assert ifc.PYRATE_ALPHA in ifg.meta_data

    def test_save_cvd_data(self):
        from os.path import join, basename, isfile
        for ifg in self.ifgs:
            if not ifg.is_open:
                ifg.open()
            data_file = join(self.params[cf.TMPDIR],
                             'cvd_data_{b}.npy'.format(b=basename(ifg.data_path).split('.')[0]))
            assert isfile(data_file)

if __name__ == "__main__":
    unittest.main()
