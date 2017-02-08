"""
Collection of tests for validating PyRate's time series analysis code.
"""
from __future__ import print_function
import os
import shutil
import sys
import tempfile
import unittest
import pytest
from datetime import date, timedelta
from numpy import nan, asarray, where
import numpy as np
from numpy.testing import assert_array_almost_equal

from pyrate import config as cf
from pyrate import mst
from pyrate import ref_phs_est as rpe
from pyrate import shared
from pyrate import vcm
from pyrate import mpiops
from pyrate import refpixel
from pyrate.config import PARALLEL, PROCESSES, NO_DATA_VALUE
from pyrate.config import TIME_SERIES_INTERP, TIME_SERIES_PTHRESH, \
    NAN_CONVERSION
from pyrate.config import TIME_SERIES_SM_FACTOR, TIME_SERIES_METHOD
from pyrate.config import TIME_SERIES_SM_ORDER
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate.timeseries import time_series
from tests.common import SYD_TEST_DIR
from tests.common import sydney_data_setup
from tests import common


def default_params():
    return {TIME_SERIES_METHOD: 1,
            TIME_SERIES_PTHRESH: 0,
            TIME_SERIES_SM_ORDER: 2,
            TIME_SERIES_SM_FACTOR: -0.25,
            TIME_SERIES_INTERP: 1,
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
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        refx, refy = run_pyrate.ref_pixel_calc(dest_paths, params)
        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)
        ifgs = shared.prepare_ifgs_without_phase(dest_paths, params)
        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        maxvar = [vcm.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 1
        params[cf.PARALLEL] = 0
        # Calculate time series
        cls.tsincr_0, cls.tscum_0, _ = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 1
        cls.tsincr_1, cls.tscum_1, cls.tsvel_1 = \
            run_pyrate.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        cls.tsincr_2, cls.tscum_2, cls.tsvel_2 = \
            run_pyrate.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

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
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)

        refx, refy = run_pyrate.ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)
        ifgs = shared.prepare_ifgs_without_phase(dest_paths, params)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        # Calculate interferogram noise
        maxvar = [vcm.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 2
        params[cf.PARALLEL] = 1
        # Calculate time series
        cls.tsincr, cls.tscum, _ = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2

        # Calculate time series
        cls.tsincr_2, cls.tscum_2, _ = \
            run_pyrate.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 0
        # Calculate time series serailly by the pixel
        cls.tsincr_0, cls.tscum_0, _ = \
            run_pyrate.calculate_time_series(ifgs, params, vcmt, mst=mst_grid)

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


@pytest.fixture(params=range(1, 6))
def modify_config(request, tempdir, get_config):
    test_conf = common.SYDNEY_TEST_CONF
    params_dict = get_config(test_conf)
    params_dict[cf.IFG_LKSX] = request.param
    params_dict[cf.IFG_LKSY] = request.param
    params_dict[cf.OBS_DIR] = tempdir()
    shared.copytree(common.SYD_TEST_GAMMA, params_dict[cf.OBS_DIR])
    params_dict[cf.IFG_FILE_LIST] = os.path.join(
        params_dict[cf.OBS_DIR], 'ifms_17')
    params_dict[cf.PARALLEL] = False
    params_dict[cf.APS_CORRECTION] = 0
    yield params_dict
    # clean up
    shutil.rmtree(params_dict[cf.OBS_DIR])


def test_timeseries_mpi(mpisync, tempdir, modify_config, ref_est_method,
                        row_splits, col_splits):
    params_dict = modify_config
    outdir = mpiops.run_once(tempdir)
    params_dict[cf.OUT_DIR] = outdir
    params_dict[cf.REF_EST_METHOD] = ref_est_method
    xlks, ylks, crop = cf.transform_params(params_dict)
    print("xlks, row_splits, col_splits, rank")
    print(xlks, row_splits, col_splits, mpiops.rank)
    if xlks * col_splits > 45 or ylks * row_splits > 70:
        print('skipping test')
        return

    base_unw_paths = cf.original_ifg_paths(params_dict[cf.IFG_FILE_LIST])
    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params_dict, xlks)

    # run prepifg, create the dest_paths files
    if mpiops.rank == 0:
        run_prepifg.gamma_prepifg(base_unw_paths, params_dict)

    mpiops.comm.barrier()

    tiles = mpiops.run_once(run_pyrate.get_tiles, dest_paths[0],
                            rows=row_splits, cols=col_splits)
    preread_ifgs = run_pyrate.create_ifg_dict(dest_paths,
                                              params=params_dict,
                                              tiles=tiles)
    run_pyrate.mpi_mst_calc(dest_paths, params_dict, tiles, preread_ifgs)
    refpx, refpy = run_pyrate.ref_pixel_calc(dest_paths, params_dict)
    run_pyrate.orb_fit_calc(dest_paths, params_dict)
    run_pyrate.ref_phase_estimation_mpi(dest_paths, params_dict, refpx, refpy)

    maxvar, vcmt = run_pyrate.maxvar_vcm_mpi(dest_paths, params_dict,
                                             preread_ifgs)
    run_pyrate.save_numpy_phase(dest_paths, tiles, params_dict)
    run_pyrate.time_series_mpi(dest_paths, params_dict, vcmt, tiles,
                               preread_ifgs)
    ifgs_mpi_out_dir = params_dict[cf.OUT_DIR]
    ifgs_mpi = sydney_data_setup(datafiles=dest_paths)

    # old vcm/maxvar estimate
    params_dict_old = modify_config
    params_dict_old[cf.OUT_DIR] = tempdir()
    params_dict_old[cf.REF_EST_METHOD] = ref_est_method
    if mpiops.rank == 0:
        xlks, ylks, crop = cf.transform_params(params_dict_old)
        base_unw_paths = cf.original_ifg_paths(
            params_dict_old[cf.IFG_FILE_LIST])
        dest_paths = cf.get_dest_paths(
            base_unw_paths, crop, params_dict_old, xlks)
        run_prepifg.gamma_prepifg(base_unw_paths, params_dict_old)

        ifgs = shared.pre_prepare_ifgs(dest_paths, params_dict_old)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params_dict_old)
        refy, refx = refpixel.ref_pixel(ifgs, params_dict)

        run_pyrate.remove_orbital_error(ifgs, params_dict)
        ifgs = shared.prepare_ifgs_without_phase(dest_paths, params_dict)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params_dict, refx, refy)
        maxvar_s = [vcm.cvd(i, params_dict_old)[0] for i in ifgs]
        vcmt_s = vcm.get_vcmt(ifgs, maxvar)
        tsincr, tscum, _ = run_pyrate.calculate_time_series(
            ifgs, params_dict_old, vcmt_s, mst=mst_grid)

        mst_mpi = reconstruct_mst(ifgs[0].shape, tiles, ifgs_mpi_out_dir)
        np.testing.assert_array_almost_equal(mst_grid, mst_mpi)
        tsincr_mpi, tscum_mpi = reconstruct_times_series(ifgs[0].shape,
                                                         tiles,
                                                         ifgs_mpi_out_dir)
        np.testing.assert_array_almost_equal(maxvar, maxvar_s)
        np.testing.assert_array_almost_equal(vcmt, vcmt_s)
        for i, j in zip(ifgs, ifgs_mpi):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data)
        np.testing.assert_array_almost_equal(tsincr, tsincr_mpi, decimal=4)
        np.testing.assert_array_almost_equal(tscum, tscum_mpi, decimal=4)
        shutil.rmtree(ifgs_mpi_out_dir)  # remove mpi out dir
        shutil.rmtree(params_dict_old[cf.OUT_DIR])  # remove serial out dir


def reconstruct_times_series(shape, tiles, output_dir):
    tsincr_file_0 = os.path.join(output_dir, 'tsincr_{}.npy'.format(0))
    shape3 = np.load(tsincr_file_0).shape[2]

    tsincr_mpi = np.empty(shape=(shape + (shape3,)), dtype=np.float32)
    tscum_mpi = np.empty_like(tsincr_mpi, dtype=np.float32)

    for i, t in enumerate(tiles):
        tsincr_file_n = os.path.join(output_dir,
                                     'tsincr_{}.npy'.format(i))
        tsincr_mpi[t.top_left_y:t.bottom_right_y,
                   t.top_left_x: t.bottom_right_x, :] = np.load(tsincr_file_n)

        tscum_file_n = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))

        tscum_mpi[t.top_left_y:t.bottom_right_y,
                  t.top_left_x: t.bottom_right_x, :] = np.load(tscum_file_n)

    return tsincr_mpi, tscum_mpi


def reconstruct_mst(shape, tiles, output_dir):
    mst_file_0 = os.path.join(output_dir, 'mst_mat_{}.npy'.format(0))
    shape0 = np.load(mst_file_0).shape[0]

    mst_mpi = np.empty(shape=((shape0,) + shape), dtype=np.float32)
    mst_mpi = np.empty_like(mst_mpi, dtype=np.float32)

    for i, t in enumerate(tiles):
        mst_file_n = os.path.join(output_dir, 'mst_mat_{}.npy'.format(i))
        mst_mpi[:, t.top_left_y:t.bottom_right_y,
                t.top_left_x: t.bottom_right_x] = np.load(mst_file_n)
    return mst_mpi


if __name__ == "__main__":
    unittest.main()
