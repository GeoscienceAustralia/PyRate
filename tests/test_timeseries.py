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
This Python module contains tests for the timeseries.py PyRate module.
"""
from __future__ import print_function
import os
import shutil
import sys
import tempfile
import unittest
from datetime import date, timedelta
from numpy import nan, asarray, where
import numpy as np
from numpy.testing import assert_array_almost_equal

import pyrate.orbital
import tests.common
from pyrate import config as cf
from pyrate import mst
from pyrate import ref_phs_est as rpe
from pyrate import shared
from pyrate import vcm
from pyrate.config import PARALLEL, PROCESSES, NO_DATA_VALUE
from pyrate.config import TIME_SERIES_PTHRESH, NAN_CONVERSION
from pyrate.config import TIME_SERIES_SM_FACTOR, TIME_SERIES_METHOD
from pyrate.config import TIME_SERIES_SM_ORDER
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate.timeseries import time_series
from tests.common import SYD_TEST_DIR, prepare_ifgs_without_phase
from tests.common import sydney_data_setup


def default_params():
    return {TIME_SERIES_METHOD: 1,
            TIME_SERIES_PTHRESH: 0,
            TIME_SERIES_SM_ORDER: 2,
            TIME_SERIES_SM_FACTOR: -0.25,
            PARALLEL: 0,
            PROCESSES: 1,
            NAN_CONVERSION: 1,
            NO_DATA_VALUE: 0}


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
        cls.ifgs = sydney_data_setup()
        cls.params = default_params()
        cls.mstmat = mst.mst_boolean_array(cls.ifgs)
        cls.maxvar = [vcm.cvd(i, cls.params)[0] for i in cls.ifgs]
        cls.vcmt = vcm.get_vcmt(cls.ifgs, cls.maxvar)

    # def test_time_series(self):
    # Sudipta:
    # 1. This has been replaced due to change in the time_series code.
    # 2. See MatlabTimeSeriesEquality for a more comprehensive test of the
    # new function.

    #     """
    #     Checks that the code works the same as the pirate MatLab code
    #     """
    #     tsincr, tscum, tsvel = time_series(
    #         self.ifgs, pthresh=self.params[TIME_SERIES_PTHRESH],
    #         params=self.params, vcmt=self.vcmt, mst=self.mstmat)
    #     expected = asarray([
    #         -11.09124207, -2.24628582, -11.37726666, -7.98105646,
    #         -8.95696049, -4.35343281, -10.64072681, -2.56493343,
    #         -6.24070214, -7.64697381, -8.59650367, -9.55619863])
    #
    #     assert_array_almost_equal(tscum[10, 10, :], expected)

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


class MatlabTimeSeriesEquality(unittest.TestCase):
    """
    Checks the python function to that of Matlab pirate ts.m and tsinvlap.m
    functionality.
    """

    @classmethod
    def setUpClass(cls):
        params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))

        cls.temp_out_dir = tempfile.mkdtemp()

        sys.argv = ['run_prepifg.py', os.path.join(SYD_TEST_DIR,
                                                   'pyrate_system_test.conf')]
        params[cf.OUT_DIR] = cls.temp_out_dir
        run_prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop, params, xlks)
        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = tests.common.mst_calculation(dest_paths, params)
        refx, refy = run_pyrate.ref_pixel_calc(dest_paths, params)
        # Estimate and remove orbit errors
        pyrate.orbital.remove_orbital_error(ifgs, params)
        ifgs = prepare_ifgs_without_phase(dest_paths, params)
        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        maxvar = [vcm.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 1
        params[cf.PARALLEL] = 0
        # Calculate time series
        cls.tsincr_0, cls.tscum_0, _ = tests.common.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 1
        cls.tsincr_1, cls.tscum_1, cls.tsvel_1 = \
            tests.common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        cls.tsincr_2, cls.tscum_2, cls.tsvel_2 = \
            tests.common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        # load the matlab data
        syd_ts_dir = os.path.join(SYD_TEST_DIR, 'matlab_time_series')
        tsincr_path = os.path.join(syd_ts_dir,
                                   'ts_incr_interp0_method1.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        # the matlab tsvel return is a bit pointless and not tested here
        # tserror is not returned
        # tserr_path = os.path.join(SYD_TIME_SERIES_DIR,
        # 'ts_error_interp0_method1.csv')
        # ts_err = np.genfromtxt(tserr_path, delimiter=',')
        tscum_path = os.path.join(syd_ts_dir,
                                  'ts_cum_interp0_method1.csv')
        ts_cum = np.genfromtxt(tscum_path)
        cls.ts_incr = np.reshape(ts_incr,
                                 newshape=cls.tsincr_0.shape, order='F')
        cls.ts_cum = np.reshape(ts_cum, newshape=cls.tscum_0.shape, order='F')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

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

    def test_time_series_equality_parallel_by_the_pixel(self):
        """
        check time series parallel by pixel jobs
        """
        self.assertEqual(self.tsincr_2.shape, self.tscum_2.shape)
        self.assertEqual(self.tsvel_2.shape, self.tsincr_2.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_2, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_2, decimal=3)

    def test_time_series_equality_serial_by_the_pixel(self):
        """
        check time series
        """

        self.assertEqual(self.tsincr_0.shape, self.tscum_0.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_0, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_0, decimal=3)


class MatlabTimeSeriesEqualityMethod2Interp0(unittest.TestCase):
    """
    Checks the python function to that of Matlab pirate ts.m and tsinvlap.m
    functionality.
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

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop, params, xlks)
        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = tests.common.mst_calculation(dest_paths, params)

        refx, refy = run_pyrate.ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.orbital.remove_orbital_error(ifgs, params)
        ifgs = prepare_ifgs_without_phase(dest_paths, params)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        # Calculate interferogram noise
        maxvar = [vcm.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 2
        params[cf.PARALLEL] = 1
        # Calculate time series
        cls.tsincr, cls.tscum, _ = tests.common.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2

        # Calculate time series
        cls.tsincr_2, cls.tscum_2, _ = \
            tests.common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 0
        # Calculate time series serailly by the pixel
        cls.tsincr_0, cls.tscum_0, _ = \
            tests.common.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        # copy matlab data
        SYD_TIME_SERIES_DIR = os.path.join(SYD_TEST_DIR, 'matlab_time_series')
        tsincr_path = os.path.join(SYD_TIME_SERIES_DIR,
                                   'ts_incr_interp0_method2.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        # the matlab tsvel return is a bit pointless and not tested here
        # tserror is not returned
        # tserr_path = os.path.join(SYD_TIME_SERIES_DIR,
        # 'ts_error_interp0_method1.csv')
        # ts_err = np.genfromtxt(tserr_path, delimiter=',')
        tscum_path = os.path.join(SYD_TIME_SERIES_DIR,
                                  'ts_cum_interp0_method2.csv')
        ts_cum = np.genfromtxt(tscum_path)

        cls.ts_incr = np.reshape(ts_incr,
                                 newshape=cls.tsincr_0.shape, order='F')
        cls.ts_cum = np.reshape(ts_cum, newshape=cls.tscum_0.shape, order='F')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_time_series_equality_parallel_by_rows(self):

        self.assertEqual(self.tsincr.shape, self.tscum.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr, decimal=1)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum, decimal=1)

    def test_time_series_equality_parallel_by_the_pixel(self):

        self.assertEqual(self.tsincr_2.shape, self.tscum_2.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_2, decimal=1)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_2, decimal=1)

    def test_time_series_equality_serial_by_the_pixel(self):

        self.assertEqual(self.tsincr_0.shape, self.tscum_0.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_0, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_0, decimal=3)

if __name__ == "__main__":
    unittest.main()
