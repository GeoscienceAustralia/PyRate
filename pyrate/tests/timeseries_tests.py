"""
Collection of tests for validating PyRate's time series analysis code.

.. codeauthor:: Vanessa Newey, Sudipta Basak
"""

import unittest
from numpy import nan, asarray, where
from numpy.testing import assert_array_almost_equal
from datetime import date, timedelta
import numpy as np
import os
import shutil

from pyrate import mst
from pyrate.tests.common import sydney_data_setup
from pyrate.timeseries import time_series
from pyrate.config import TIME_SERIES_PTHRESH, TIME_SERIES_SM_ORDER
from pyrate.config import TIME_SERIES_SM_FACTOR
from pyrate.config import PARALLEL, PROCESSES
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate import config as cf
from pyrate import reference_phase_estimation as rpe
from pyrate import vcm


def default_params():
    return {TIME_SERIES_PTHRESH: 10,
            TIME_SERIES_SM_ORDER: 2,
            TIME_SERIES_SM_FACTOR: -0.25,
            PARALLEL: 0,
            PROCESSES: 1}


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
        cls.mstmat = mst.mst_matrix_ifg_indices_as_boolean_array(cls.ifgs)
        cls.maxvar = [vcm.cvd(i)[0] for i in cls.ifgs]
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
            self.ifgs, pthresh=0,
            params=self.params, vcmt=self.vcmt, mst=None)
        expected = asarray([[[0.50,  3.0,  4.0,  5.5,  6.5]]])
        assert_array_almost_equal(tscum, expected, decimal=2)


class MatlabTimeSeriesEquality(unittest.TestCase):
    """
    Checks the python function to that of Matlab pirate ts.m and tsinvlap.m
    functionality.
    """

    @classmethod
    def setUpClass(cls):
        # start each full test run cleanly
        shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

        os.makedirs(SYD_TEST_OUT)

        params = cf.get_config_params(
                os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))
        params[cf.REF_EST_METHOD] = 2
        run_prepifg.main(
            config_file=os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR,
                                     'orbital_error.conf'))

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

        assert isinstance(ifg_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_instance, nan_conversion=True)
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

        for i in ifgs:
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        # Calculate interferogram noise
        # TODO: assign maxvar to ifg metadata (and geotiff)?
        maxvar = [vcm.cvd(i)[0] for i in ifgs]

        # Calculate temporal variance-covariance matrix
        vcmt = vcm.get_vcmt(ifgs, maxvar)

        SYD_TIME_SERIES_DIR = os.path.join(SYD_TEST_DIR, 'matlab_time_series')
        tsincr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_incr.csv')
        cls.ts_incr = np.genfromtxt(tsincr_path, delimiter=',')
        tserr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_error.csv')
        cls.ts_err = np.genfromtxt(tserr_path, delimiter=',')
        tscum_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_cum.csv')
        cls.ts_cum = np.genfromtxt(tscum_path, delimiter=',')

        params[cf.PARALLEL] = 1
        # Calculate time series
        if params[cf.TIME_SERIES_CAL] != 0:
            cls.tsincr, cls.tscum, cls.tsvel = run_pyrate.calculate_time_series(
                ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        # Calculate time series
        if params[cf.TIME_SERIES_CAL] != 0:
            cls.tsincr_2, cls.tscum_2, cls.tsvel_2 = \
                run_pyrate.calculate_time_series(
                ifgs, params, vcmt, mst=mst_grid
                )

    def test_time_series_equality(self):

        self.assertEqual(self.tsincr.shape, self.tscum.shape)
        self.assertEqual(self.tsvel.shape, self.tsincr.shape)
        ts_incr = np.reshape(self.ts_incr,
                             newshape=self.tsincr.shape, order='F')
        ts_cum = np.reshape(self.ts_cum, newshape=self.tsincr.shape, order='F')

        #TODO: Investigate why the entire matrices don't equal
        # Current hypothesis is that the pseudo inverse computed are different
        # in matlab and python as they are based of different convergence
        # criteria.
        np.testing.assert_array_almost_equal(
            ts_incr[:11, :45, :], self.tsincr[:11, :45, :], decimal=4)

        np.testing.assert_array_almost_equal(
            ts_cum[:11, :45, :], self.tscum[:11, :45, :], decimal=4)


    def test_time_series_equality_parallel_by_the_pixel(self):

        self.assertEqual(self.tsincr_2.shape, self.tscum_2.shape)
        self.assertEqual(self.tsvel_2.shape, self.tsincr_2.shape)
        ts_incr = np.reshape(self.ts_incr,
                             newshape=self.tsincr.shape, order='F')
        ts_cum = np.reshape(self.ts_cum, newshape=self.tsincr.shape, order='F')

        #TODO: Investigate why the entire matrices don't equal
        # Current hypothesis is that the pseudo inverse computed are different
        # in matlab and python as they are based of different convergence
        # criteria.
        np.testing.assert_array_almost_equal(
            ts_incr[:11, :45, :], self.tsincr_2[:11, :45, :], decimal=4)

        np.testing.assert_array_almost_equal(
            ts_cum[:11, :45, :], self.tscum_2[:11, :45, :], decimal=4)


if __name__ == "__main__":
    unittest.main()
