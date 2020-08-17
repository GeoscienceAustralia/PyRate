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
This Python module contains tests for the timeseries.py PyRate module.
"""
import os
import shutil
import pytest
from datetime import date, timedelta
from numpy import nan, asarray, where, array
import numpy as np
from numpy.testing import assert_array_almost_equal

import pyrate.core.orbital
import pyrate.core.ref_phs_est
import pyrate.core.refpixel
import tests.common as common
from pyrate.core import config as cf, mst, covariance
from pyrate import correct, prepifg, conv2tif
from pyrate.configuration import Configuration
from pyrate.core.timeseries import time_series, linear_rate_pixel, linear_rate_array, TimeSeriesError


def default_params():
    return {cf.TIME_SERIES_METHOD: 1,
            cf.TIME_SERIES_PTHRESH: 0,
            cf.TIME_SERIES_SM_ORDER: 2,
            cf.TIME_SERIES_SM_FACTOR: -0.25,
            cf.PARALLEL: 0,
            cf.PROCESSES: 1,
            cf.NAN_CONVERSION: 1,
            cf.NO_DATA_VALUE: 0}


class SinglePixelIfg(object):
    """
    A single pixel ifg (interferogram) solely for unit testing
    """

    def __init__(self, first, second, phase, nan_fraction):
        self.phase_data = asarray([[phase]])
        self.first = first
        self.second = second
        self.nrows = 1
        self.ncols = 1
        self.nan_fraction = asarray([nan_fraction])

    def convert_to_nans(self, val=0):
        """
        Converts given values in phase data to NaNs
        val - value to convert, default is 0
        """
        self.phase_data = where(self.phase_data == val, nan, self.phase_data)
        self.nan_converted = True


class TestTimeSeries:
    """Verifies error checking capabilities of the time_series function"""

    @classmethod
    def setup_class(cls):
        cls.ifgs = common.small_data_setup()
        cls.params = default_params()
        cls.mstmat = mst.mst_boolean_array(cls.ifgs)
        r_dist = covariance.RDist(cls.ifgs[0])()
        cls.maxvar = [covariance.cvd(i, cls.params, r_dist)[0]
                      for i in cls.ifgs]
        cls.vcmt = covariance.get_vcmt(cls.ifgs, cls.maxvar)

    def test_time_series_unit(self):
        """
        Checks that the code works the same as the calculated example
        """
        ifirst = asarray([1, 1, 2, 2, 3, 3, 4, 5])
        isecond = asarray([2, 4, 3, 4, 5, 6, 6, 6])
        timeseries = asarray([0.0, 0.1, 0.6, 0.8, 1.1, 1.3])
        phase = asarray([0.5, 4, 2.5, 3.5, 2.5, 3.5, 2.5, 1])
        nan_fraction = asarray([0.5, 0.4, 0.2, 0.3, 0.1, 0.3, 0.2, 0.1])

        now = date.today()

        dates = [now + timedelta(days=(t*365.25)) for t in timeseries]
        dates.sort()
        first = [dates[m_num - 1] for m_num in ifirst]
        second = [dates[s_num - 1] for s_num in isecond]

        self.ifgs = [SinglePixelIfg(m, s, p, n) for m, s, p, n in
                     zip(first, second, phase, nan_fraction)]

        tsincr, tscum, tsvel = time_series(
            self.ifgs, params=self.params, vcmt=self.vcmt, mst=None)
        expected = asarray([[[0.50, 3.0, 4.0, 5.5, 6.5]]])
        assert_array_almost_equal(tscum, expected, decimal=2)


class TestLegacyTimeSeriesEquality:

    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        params[cf.TEMP_MLOOKED_DIR] = os.path.join(params[cf.OUT_DIR], cf.TEMP_MLOOKED_DIR)
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, _, crop = cf.transform_params(params)

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)
        correct._copy_mlooked(params)
        copied_dest_paths = [os.path.join(params[cf.TEMP_MLOOKED_DIR], os.path.basename(d)) for d in dest_paths]
        del dest_paths
        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(copied_dest_paths, params)
        mst_grid = common.mst_calculation(copied_dest_paths, params)
        refx, refy = pyrate.core.refpixel.ref_pixel_calc_wrapper(params)

        params[cf.REFX] = refx
        params[cf.REFY] = refy
        params[cf.ORBFIT_OFFSET] = True
        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)
        ifgs = common.prepare_ifgs_without_phase(copied_dest_paths, params)
        for ifg in ifgs:
            ifg.close()
        correct._update_params_with_tiles(params)
        _, ifgs = pyrate.core.ref_phs_est.ref_phase_est_wrapper(params)
        ifgs[0].open()
        r_dist = covariance.RDist(ifgs[0])()
        ifgs[0].close()
        maxvar = [covariance.cvd(i, params, r_dist)[0] for i in copied_dest_paths]
        for ifg in ifgs:
            ifg.open()
        vcmt = covariance.get_vcmt(ifgs, maxvar)

        for ifg in ifgs:
            ifg.close()
            ifg.open()
            ifg.nodata_value = 0.0

        params[cf.TIME_SERIES_METHOD] = 1
        params[cf.PARALLEL] = 0
        # Calculate time series
        cls.tsincr_0, cls.tscum_0, _ = common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 1
        cls.tsincr_1, cls.tscum_1, cls.tsvel_1 = common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        # load the legacy data
        ts_dir = os.path.join(common.SML_TEST_DIR, 'time_series')
        tsincr_path = os.path.join(ts_dir, 'ts_incr_interp0_method1.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        tscum_path = os.path.join(ts_dir, 'ts_cum_interp0_method1.csv')
        ts_cum = np.genfromtxt(tscum_path)
        cls.ts_incr = np.reshape(ts_incr, newshape=cls.tsincr_0.shape, order='F')
        cls.ts_cum = np.reshape(ts_cum, newshape=cls.tscum_0.shape, order='F')
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_time_series_equality_parallel_by_rows(self):
        """
        check time series parallel by rows jobs
        """

        self.assertEqual(self.tsincr_1.shape, self.tscum_1.shape)
        self.assertEqual(self.tsvel_1.shape, self.tsincr_1.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_1, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_1, decimal=3)

    def test_time_series_equality_serial_by_the_pixel(self):
        """
        check time series
        """

        self.assertEqual(self.tsincr_0.shape, self.tscum_0.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_0, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_0, decimal=3)

    @staticmethod
    def assertEqual(val1, val2):
        assert val1 == val2


class TestLegacyTimeSeriesEqualityMethod2Interp0:

    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        params[cf.TEMP_MLOOKED_DIR] = os.path.join(params[cf.OUT_DIR], cf.TEMP_MLOOKED_DIR)
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, _, crop = cf.transform_params(params)

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)
        correct._copy_mlooked(params)
        copied_dest_paths = [os.path.join(params[cf.TEMP_MLOOKED_DIR], os.path.basename(d)) for d in dest_paths]
        del dest_paths
        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(copied_dest_paths, params)
        mst_grid = common.mst_calculation(copied_dest_paths, params)
        refx, refy = pyrate.core.refpixel.ref_pixel_calc_wrapper(params)

        params[cf.REFX] = refx
        params[cf.REFY] = refy
        params[cf.ORBFIT_OFFSET] = True

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)
        ifgs = common.prepare_ifgs_without_phase(copied_dest_paths, params)
        for ifg in ifgs:
            ifg.close()

        correct._update_params_with_tiles(params)
        _, ifgs = pyrate.core.ref_phs_est.ref_phase_est_wrapper(params)
        ifgs[0].open()
        r_dist = covariance.RDist(ifgs[0])()
        ifgs[0].close()
        # Calculate interferogram noise
        maxvar = [covariance.cvd(i, params, r_dist)[0] for i in copied_dest_paths]
        for ifg in ifgs:
            ifg.open()
        vcmt = covariance.get_vcmt(ifgs, maxvar)
        for ifg in ifgs:
            ifg.close()
            ifg.open()
            ifg.nodata_value = 0.0

        params[cf.TIME_SERIES_METHOD] = 2
        params[cf.PARALLEL] = 1
        # Calculate time series
        cls.tsincr, cls.tscum, _ = common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 0
        # Calculate time series serailly by the pixel
        cls.tsincr_0, cls.tscum_0, _ = common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        # copy legacy data
        SML_TIME_SERIES_DIR = os.path.join(common.SML_TEST_DIR, 'time_series')
        tsincr_path = os.path.join(SML_TIME_SERIES_DIR, 'ts_incr_interp0_method2.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        tscum_path = os.path.join(SML_TIME_SERIES_DIR, 'ts_cum_interp0_method2.csv')
        ts_cum = np.genfromtxt(tscum_path)

        cls.ts_incr = np.reshape(ts_incr, newshape=cls.tsincr_0.shape, order='F')
        cls.ts_cum = np.reshape(ts_cum, newshape=cls.tscum_0.shape, order='F')
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_time_series_equality_parallel_by_rows(self):

        assert self.tsincr.shape == self.tscum.shape

        np.testing.assert_array_almost_equal(self.ts_incr, self.tsincr, decimal=1)

        np.testing.assert_array_almost_equal(self.ts_cum, self.tscum, decimal=1)

    def test_time_series_equality_serial_by_the_pixel(self):

        assert self.tsincr_0.shape == self.tscum_0.shape

        np.testing.assert_array_almost_equal(self.ts_incr, self.tsincr_0, decimal=3)

        np.testing.assert_array_almost_equal(self.ts_cum, self.tscum_0, decimal=3)


class TestLinearRatePixel:
    """
    Tests the linear regression algorithm for determining the best
    fitting velocity from a cumulative time series
    """

    def test_linear_rate_pixel_clean(self):
        y = array([0, 2, 4, 6, 8, 10])
        t = array([0, 1, 2, 3, 4, 5])
        exp = (2.0, 0.0, 1.0, 0.0, 6)
        res = linear_rate_pixel(y, t)
        assert res == exp

    def test_linear_rate_pixel_neg_rate(self):
        y = array([0, -2, -4, -6, -8, -10])
        t = array([0, 1, 2, 3, 4, 5])
        exp = (-2.0, 0.0, 1.0, 0.0, 6)
        res = linear_rate_pixel(y, t)
        assert res == exp

    def test_linear_rate_pixel_outlier(self):
        y = array([0, 2, 4, 6, 8, 20])
        t = array([0, 1, 2, 3, 4, 5])
        exp = (3.428571, -1.904761, 0.812030, 0.824786, 6)
        res = linear_rate_pixel(y, t)
        assert res == pytest.approx(exp, rel=1e-6)

    def test_linear_rate_pixel_noise(self):
        y = array([0, 2, 4, 6, 8, 10])
        r = y + np.random.rand(6) # add different uniform noise each time
        t = array([0, 1, 2, 3, 4, 5])
        exprate = 2.0
        explsqd = 1.0
        experr = 0.0
        rate, _, lsqd, err, _ = linear_rate_pixel(y, t)
        assert exprate == pytest.approx(rate, rel=1e-1)
        assert explsqd == pytest.approx(lsqd, rel=1e-1)
        assert experr == pytest.approx(err, rel=1e-1)

    def test_linear_rate_pixel_exception(self):
        # input vectors should be equal length
        y = array([2, 4, 6, 8, 10])
        t = array([0, 1, 2, 3, 4, 5])
        with pytest.raises(TimeSeriesError):
            res = linear_rate_pixel(y, t)

    def test_linear_rate_pixel_nans(self):
        # at least two obs are required for line fitting
        y = array([0, nan, nan, nan, nan, nan])
        t = array([0, 1, 2, 3, 4, 5])
        exp = (nan, nan, nan, nan, nan)
        res = linear_rate_pixel(y, t)
        assert res == exp


class TestLinearRateArray:
    """
    Tests the array loop wrapper for the linear regression algorithm using real data
    """

    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls, roipac_params):
        cls.params = roipac_params
        cls.ifgs = common.small_data_setup()

        # read in input (tscuml) and expected output arrays
        tscuml_path = os.path.join(common.SML_TEST_LINRATE, "tscuml_0.npy")
        cls.tscuml0 = np.load(tscuml_path)
        # add zero epoch to tscuml 3D array
        cls.tscuml = np.insert(cls.tscuml0, 0, 0, axis=2) 

        linrate_path = os.path.join(common.SML_TEST_LINRATE, "linear_rate.npy")
        cls.linrate = np.load(linrate_path)

        error_path = os.path.join(common.SML_TEST_LINRATE, "linear_error.npy")
        cls.error = np.load(error_path)

        icpt_path = os.path.join(common.SML_TEST_LINRATE, "linear_intercept.npy")
        cls.icpt = np.load(icpt_path)

        samp_path = os.path.join(common.SML_TEST_LINRATE, "linear_samples.npy")
        cls.samp = np.load(samp_path)

        rsq_path = os.path.join(common.SML_TEST_LINRATE, "linear_rsquared.npy")
        cls.rsq = np.load(rsq_path)

    def test_linear_rate_array(self):
        """
        Input and expected output are on disk. This test only tests the linear_rate_array
        and linear_rate_pixel functions using real data.
        """
        l, i, r, e, s = linear_rate_array(self.tscuml, self.ifgs, self.params)
        # test to 20 decimal places
        assert_array_almost_equal(self.linrate, l, 1e-20)
        assert_array_almost_equal(self.icpt, i, 1e-20)
        assert_array_almost_equal(self.rsq, r, 1e-20)
        assert_array_almost_equal(self.error, e, 1e-20)
        assert_array_almost_equal(self.samp, s, 1e-20)

    def test_linear_rate_array_exception(self):
        # depth of tscuml should equal nepochs
        with pytest.raises(TimeSeriesError):
            res = linear_rate_array(self.tscuml0, self.ifgs, self.params)

