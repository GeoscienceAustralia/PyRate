#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
This Python module contains tests for the APS PyRate module.
"""
import os
import shutil
import pytest
import numpy as np
from os.path import join
from scipy.ndimage import gaussian_filter1d, gaussian_filter
from numpy.testing import assert_array_almost_equal

import pyrate.constants as C
from pyrate import conv2tif, prepifg, correct
from pyrate.configuration import Configuration, MultiplePaths
from pyrate.core import shared
from pyrate.core.aps import spatio_temporal_filter, _interpolate_nans_2d, _kernel
from pyrate.core.aps import gaussian_temporal_filter as tlpfilter, gaussian_spatial_filter as slpfilter
from pyrate.core.shared import Ifg
from pyrate.core.ifgconstants import DAYS_PER_YEAR
from tests import common
from tests.common import BASE_TEST, PY37GDAL304


@pytest.fixture(params=["linear", "nearest", "cubic"])
def slpnanfill_method(request):
    return request.param


@pytest.fixture(params=[0.1, 0.5, 1, 5])
def slpfcutoff_method(request):
    return request.param


@pytest.fixture(params=[0.001, 0.01])
def slpfcutoff(request):
    return request.param


def test_interpolate_nans_2d(slpnanfill_method):
    arr = np.random.rand(200, 100)
    arr[arr < 0.1] = np.nan  # insert some nans
    assert np.sum(np.isnan(arr)) != 0  # some nans present
    _interpolate_nans_2d(arr, method=slpnanfill_method)
    assert np.sum(np.isnan(arr)) == 0  # should not be any nans


class TestSpatialFilter:
    """
    Test the implementation of Gaussian spatial filter
    """

    def setup_method(self):
        ifg_path = join(str(BASE_TEST), 'cropB', '20180106-20180130_ifg.tif')
        ifg = Ifg(ifg_path)
        ifg.open()
        p = ifg.phase_data
        self.x_size = ifg.x_size
        self.y_size = ifg.y_size
        # convert zeros to NaNs
        p[p == 0] = np.nan
        self.phase = p

    def test_gaussian_filter(self, slpfcutoff_method):
        """
        Compare against scipy.ndimage.gaussian_filter,
        that operates in the spatial domain on data with equal resolution in x and y
        """
        e = gaussian_filter(self.phase, sigma=slpfcutoff_method)
        r = slpfilter(self.phase, cutoff=slpfcutoff_method,
                      x_size=self.x_size, y_size=self.y_size)
        # pull data from middle of images; away from potential edge effects
        exp = e[50:150, 50:150]
        res = r[50:150, 50:150]
        assert_array_almost_equal(res, exp, 0)


def test_gaussian_kernel():
    """
    Test the Gaussian smoothing kernel
    """
    x = np.arange(1, 10, 2)
    res = _kernel(x, 3)
    exp = np.array([0.94595947, 0.60653066, 0.24935221, 0.06572853, 0.011109])
    np.testing.assert_array_almost_equal(res, exp, decimal=6)

    res = _kernel(x, 4)
    exp = np.array([0.96923323, 0.7548396, 0.45783336, 0.21626517, 0.07955951])
    np.testing.assert_array_almost_equal(res, exp, decimal=6)

    res = _kernel(x, 0.001)
    exp = np.array([0, 0, 0, 0, 0])
    np.testing.assert_array_equal(res, exp)


class TestTemporalFilter:
    """
    Tests for the temporal filter with synthetic data for a single pixel
    """

    def setup_method(self):
        self.thr = 1  # no nans in these test cases, threshold = 1
        # instance of normally distributed noise
        n = np.array([-0.36427456, 0.69539061, 0.42181139, -2.56306134,
                      0.55844095, -0.65562626, 0.65607911, 1.19431637,
                      -1.43837395, -0.91656358])
        # synthetic incremental displacement
        d = np.array([1., 1., 0.7, 0.3, 0., 0.1, 0.2, 0.6, 1., 1.])
        # incremental displacement + noise
        self.tsincr = d * 2 + n
        # regular time series, every 12 days
        self.interval = 12 / DAYS_PER_YEAR  # 0.03285 years
        intv = np.ones(d.shape, dtype=np.float32) * self.interval
        self.span = np.cumsum(intv)

    def test_tlpfilter_repeatability(self):
        """
        TEST 1: check for repeatability against expected result from
        tlpfilter. Cutoff equal to sampling interval (sigma=1)
        """
        res = tlpfilter(self.tsincr, self.interval, self.span, self.thr)
        exp = np.array([1.9936507, 1.9208364, 1.0252733, -0.07402889,
                        -0.1842336, 0.24325351, 0.94737214, 1.3890865,
                        1.1903466, 1.0036403])
        np.testing.assert_array_almost_equal(res, exp, decimal=6)

    def test_tlpfilter_scipy_sig1(self):
        """
        TEST 2: compare tlpfilter to scipy.ndimage.gaussian_filter1d. Works for
        regularly sampled data. Cutoff equal to sampling interval (sigma=1)
        """
        res = tlpfilter(self.tsincr, self.interval, self.span, self.thr)
        exp = gaussian_filter1d(self.tsincr, sigma=1)
        np.testing.assert_array_almost_equal(res, exp, decimal=1)

    def test_tlpfilter_scipy_sig2(self):
        """
        TEST 3: compare tlpfilter to scipy.ndimage.gaussian_filter1d. Works for
        regularly sampled data. Cutoff equal to twice the sampling interval (sigma=2)
        """
        res = tlpfilter(self.tsincr, self.interval * 2, self.span, self.thr)
        exp = gaussian_filter1d(self.tsincr, sigma=2)
        np.testing.assert_array_almost_equal(res, exp, decimal=1)

    def test_tlpfilter_scipy_sig05(self):
        """
        TEST 4: compare tlpfilter to scipy.ndimage.gaussian_filter1d. Works for
        regularly sampled data. Cutoff equal to half the sampling interval (sigma=0.5)
        """
        res = tlpfilter(self.tsincr, self.interval * 0.5, self.span, self.thr)
        exp = gaussian_filter1d(self.tsincr, sigma=0.5)
        np.testing.assert_array_almost_equal(res, exp, decimal=2)


@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL304, reason="Only run in one CI env")
class TestAPSErrorCorrectionsOnDiscReused:

    @classmethod
    def setup_method(cls):
        cls.conf = common.TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._update_params_with_tiles(cls.params)
        correct._create_ifg_dict(cls.params)
        multi_paths = cls.params[C.INTERFEROGRAM_FILES]
        cls.ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        cls.ifgs = [shared.Ifg(i) for i in cls.ifg_paths]
        for i in cls.ifgs:
            i.open()
        shared.save_numpy_phase(cls.ifg_paths, cls.params)
        correct.mst_calc_wrapper(cls.params)

    @classmethod
    def teardown_method(cls):
        shutil.rmtree(cls.params[C.OUT_DIR])

    def test_aps_error_files_on_disc(self, slpnanfill_method, slpfcutoff):
        self.params[C.SLPF_NANFILL_METHOD] = slpnanfill_method
        self.params[C.SLPF_CUTOFF] = slpfcutoff
        spatio_temporal_filter(self.params)

        # test_orb_errors_written
        aps_error_files = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        assert all(p.exists() for p in aps_error_files)
        last_mod_times = [os.stat(o).st_mtime for o in aps_error_files]
        phase_prev = [i.phase_data for i in self.ifgs]

        # run aps error removal again
        spatio_temporal_filter(self.params)
        aps_error_files2 = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        # if files are written again - times will change
        last_mod_times_2 = [os.stat(o).st_mtime for o in aps_error_files2]
        # test_aps_error_reused_if_params_unchanged
        assert all(a == b for a, b in zip(last_mod_times, last_mod_times_2))
        phase_now = [i.phase_data for i in self.ifgs]

        # run aps error correction once mroe
        spatio_temporal_filter(self.params)
        aps_error_files3 = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        last_mod_times_3 = [os.stat(o).st_mtime for o in aps_error_files3]
        assert all(a == b for a, b in zip(last_mod_times, last_mod_times_3))
        phase_again = [i.phase_data for i in self.ifgs]
        np.testing.assert_array_equal(phase_prev, phase_now)
        np.testing.assert_array_equal(phase_prev, phase_again)
