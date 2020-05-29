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
import sys
import tempfile
import unittest
from pathlib import Path
from datetime import date, timedelta
from numpy import nan, asarray, where
import numpy as np
from numpy.testing import assert_array_almost_equal

import pyrate.core.orbital
import tests.common as common
from pyrate.core import ref_phs_est as rpe, config as cf, mst, covariance, roipac
from pyrate import process, prepifg, conv2tif
from pyrate.configuration import Configuration
from pyrate.core.timeseries import time_series


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

    def __init__(self, master, slave, phase, nan_fraction):
        self.phase_data = asarray([[phase]])
        self.master = master
        self.slave = slave
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


class TimeSeriesTests(unittest.TestCase):
    """Verifies error checking capabilities of the time_series function"""

    @classmethod
    def setUpClass(cls):
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
        imaster = asarray([1, 1, 2, 2, 3, 3, 4, 5])
        islave = asarray([2, 4, 3, 4, 5, 6, 6, 6])
        timeseries = asarray([0.0, 0.1, 0.6, 0.8, 1.1, 1.3])
        phase = asarray([0.5, 4, 2.5, 3.5, 2.5, 3.5, 2.5, 1])
        nan_fraction = asarray([0.5, 0.4, 0.2, 0.3, 0.1, 0.3, 0.2, 0.1])

        now = date.today()

        dates = [now + timedelta(days=(t*365.25)) for t in timeseries]
        dates.sort()
        master = [dates[m_num - 1] for m_num in imaster]
        slave = [dates[s_num - 1] for s_num in islave]

        self.ifgs = [SinglePixelIfg(m, s, p, n) for m, s, p, n in
                     zip(master, slave, phase, nan_fraction)]

        tsincr, tscum, tsvel = time_series(
            self.ifgs, params=self.params, vcmt=self.vcmt, mst=None)
        expected = asarray([[[0.50, 3.0, 4.0, 5.5, 6.5]]])
        assert_array_almost_equal(tscum, expected, decimal=2)


class LegacyTimeSeriesEquality(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', common.TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        base_ifg_paths = [c.unwrapped_path for c in params[cf.INTERFEROGRAM_FILES]]
        headers = [roipac.roipac_header(i, params) for i in base_ifg_paths]
        dest_paths = [Path(cls.temp_out_dir).joinpath(Path(c.sampled_path).name).as_posix()
                      for c in params[cf.INTERFEROGRAM_FILES][:-2]]
        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        refx, refy = process._ref_pixel_calc(dest_paths, params)
        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params, headers)
        ifgs = common.prepare_ifgs_without_phase(dest_paths, params)
        for ifg in ifgs:
            ifg.close()
        _, ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)
        ifgs[0].open()
        r_dist = covariance.RDist(ifgs[0])()
        ifgs[0].close()
        maxvar = [covariance.cvd(i, params, r_dist)[0] for i in dest_paths]
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

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)
        common.remove_tifs(
            cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])

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


class LegacyTimeSeriesEqualityMethod2Interp0(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.PARALLEL] = 0
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        base_ifg_paths = [c.unwrapped_path for c in params[cf.INTERFEROGRAM_FILES]]
        headers = [roipac.roipac_header(i, params) for i in base_ifg_paths]
        dest_paths = [Path(cls.temp_out_dir).joinpath(Path(c.sampled_path).name).as_posix()
                      for c in params[cf.INTERFEROGRAM_FILES][:-2]]
        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)

        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params, headers)
        ifgs = common.prepare_ifgs_without_phase(dest_paths, params)
        for ifg in ifgs:
            ifg.close()
        _, ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)
        ifgs[0].open()
        r_dist = covariance.RDist(ifgs[0])()
        ifgs[0].close()
        # Calculate interferogram noise
        maxvar = [covariance.cvd(i, params, r_dist)[0] for i in dest_paths]
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

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.temp_out_dir)
        common.remove_tifs(cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])

    def test_time_series_equality_parallel_by_rows(self):

        assert self.tsincr.shape == self.tscum.shape

        np.testing.assert_array_almost_equal(self.ts_incr, self.tsincr, decimal=1)

        np.testing.assert_array_almost_equal(self.ts_cum, self.tscum, decimal=1)

    def test_time_series_equality_serial_by_the_pixel(self):

        self.assertEqual(self.tsincr_0.shape, self.tscum_0.shape)

        np.testing.assert_array_almost_equal(self.ts_incr, self.tsincr_0, decimal=3)

        np.testing.assert_array_almost_equal(self.ts_cum, self.tscum_0, decimal=3)


if __name__ == "__main__":
    unittest.main()
