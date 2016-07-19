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
import sys
import shutil
import tempfile
import subprocess
import glob
from itertools import product
from osgeo import gdal

from pyrate import mst
from pyrate.tests.common import sydney_data_setup
from pyrate.timeseries import time_series
from pyrate.config import TIME_SERIES_PTHRESH, TIME_SERIES_SM_ORDER
from pyrate.config import TIME_SERIES_SM_FACTOR, TIME_SERIES_METHOD
from pyrate.config import TIME_SERIES_INTERP, TIME_SERIES_PTHRESH, NAN_CONVERSION
from pyrate.config import PARALLEL, PROCESSES, NO_DATA_VALUE
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate.tests.common import SYD_TEST_DIR
from pyrate import config as cf
from pyrate import reference_phase_estimation as rpe
from pyrate import vcm
from pyrate.tests import common
from pyrate import vcm as vcm_module
from pyrate import shared
from pyrate.shared import Ifg

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
        expected = asarray([[[0.50,  3.0,  4.0,  5.5,  6.5]]])
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

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params,
                                               xlks)
        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)
        ifgs = shared.prepare_ifgs_without_phase(dest_paths, params)
        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        maxvar = [vcm_module.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm_module.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 1
        params[cf.PARALLEL] = 0
        # Calculate time series
        cls.tsincr_0, cls.tscum_0, _ = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 1
        cls.tsincr_1, cls.tscum_1, cls.tsvel_1 = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        cls.tsincr_2, cls.tscum_2, cls.tsvel_2 = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        # load the matlab data
        SYD_TIME_SERIES_DIR = os.path.join(SYD_TEST_DIR, 'matlab_time_series')
        tsincr_path = os.path.join(SYD_TIME_SERIES_DIR,
                                   'ts_incr_interp0_method1.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        # the matlab tsvel return is a bit pointless and not tested here
        # tserror is not returned
        # tserr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_error_interp0_method1.csv')
        # ts_err = np.genfromtxt(tserr_path, delimiter=',')
        tscum_path = os.path.join(SYD_TIME_SERIES_DIR,
                                  'ts_cum_interp0_method1.csv')
        ts_cum = np.genfromtxt(tscum_path)
        cls.ts_incr = np.reshape(ts_incr, newshape=cls.tsincr_0.shape, order='F')
        cls.ts_cum = np.reshape(ts_cum, newshape=cls.tscum_0.shape, order='F')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_time_series_equality_parallel_by_rows(self):

        self.assertEqual(self.tsincr_1.shape, self.tscum_1.shape)
        self.assertEqual(self.tsvel_1.shape, self.tsincr_1.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_1, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_1, decimal=3)

    def test_time_series_equality_parallel_by_the_pixel(self):

        self.assertEqual(self.tsincr_2.shape, self.tscum_2.shape)
        self.assertEqual(self.tsvel_2.shape, self.tsincr_2.shape)

        np.testing.assert_array_almost_equal(
            self.ts_incr, self.tsincr_2, decimal=3)

        np.testing.assert_array_almost_equal(
            self.ts_cum, self.tscum_2, decimal=3)

    def test_time_series_equality_serial_by_the_pixel(self):

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
        params = cf.get_config_params(
                os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))

        cls.temp_out_dir = tempfile.mkdtemp()

        sys.argv = ['run_prepifg.py', os.path.join(SYD_TEST_DIR,
                                     'pyrate_system_test.conf')]
        params[cf.OUT_DIR] = cls.temp_out_dir
        run_prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params,
                                               xlks)
        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)
        ifgs = shared.prepare_ifgs_without_phase(dest_paths, params)

        _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

        # Calculate interferogram noise
        maxvar = [vcm_module.cvd(i, params)[0] for i in ifgs]
        vcmt = vcm_module.get_vcmt(ifgs, maxvar)

        params[cf.TIME_SERIES_METHOD] = 2
        params[cf.PARALLEL] = 1
        # Calculate time series
        cls.tsincr, cls.tscum, _ = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2

        # Calculate time series
        cls.tsincr_2, cls.tscum_2, _ = \
            run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid
            )

        params[cf.PARALLEL] = 0
        # Calculate time series serailly by the pixel
        cls.tsincr_0, cls.tscum_0, _ = \
            run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid
            )

        # copy matlab data
        SYD_TIME_SERIES_DIR = os.path.join(SYD_TEST_DIR, 'matlab_time_series')
        tsincr_path = os.path.join(SYD_TIME_SERIES_DIR,
                                   'ts_incr_interp0_method2.csv')
        ts_incr = np.genfromtxt(tsincr_path)

        # the matlab tsvel return is a bit pointless and not tested here
        # tserror is not returned
        # tserr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_error_interp0_method1.csv')
        # ts_err = np.genfromtxt(tserr_path, delimiter=',')
        tscum_path = os.path.join(SYD_TIME_SERIES_DIR,
                                  'ts_cum_interp0_method2.csv')
        ts_cum = np.genfromtxt(tscum_path)

        cls.ts_incr = np.reshape(ts_incr, newshape=cls.tsincr_0.shape, order='F')
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


class MPITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF
        cls.temp_dir = tempfile.mkdtemp()
        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        cls.params[cf.IFG_FILE_LIST] = os.path.join(
            common.SYD_TEST_GAMMA, 'ifms_17')
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.APS_CORRECTION] = 0
        cls.params[cf.REF_EST_METHOD] = 1
        cls.params[cf.ORBITAL_FIT_METHOD] = 1
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        cls.base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])

    @classmethod
    def process(cls):
        cls.params[cf.OUT_DIR] = cls.tif_dir
        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            cls.base_unw_paths, crop, cls.params, xlks)

        # create the dest_paths files
        run_prepifg.gamma_prepifg(cls.base_unw_paths, cls.params)

        # now create test ifgs
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

        # give the log file any name
        cls.log_file = os.path.join(cls.tif_dir, 'timeseries_mpi.log')

        # create the conf file in out_dir
        cls.conf_file = tempfile.mktemp(suffix='.conf', dir=cls.tif_dir)
        cf.write_config_file(cls.params, cls.conf_file)

        # conf file crated?
        assert os.path.exists(cls.conf_file)

        # Calc time series using MPI
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar.py ' + \
              cls.conf_file
        cmd = str.split()
        subprocess.check_call(cmd)
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar_2.py ' + \
              cls.conf_file
        cmd = str.split()
        subprocess.check_call(cmd)
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar_3.py ' + \
              cls.conf_file
        cmd = str.split()
        subprocess.check_call(cmd)

        str = 'mpirun -np 4 python pyrate/nci/postprocessing.py ' + \
              cls.conf_file
        cmd = str.split()
        subprocess.check_call(cmd)

        # sydney test data known to have a span/nvelpar of 12
        cls.tsincr_mpi = np.empty(shape=(cls.ifgs[0].shape + (12,)),
                           dtype=np.float32)
        cls.tscum_mpi = np.empty_like(cls.tsincr_mpi, dtype=np.float32)

        tiles = shared.create_tiles(cls.ifgs[0].shape,
                                    n_ifgs=17, nrows=3, ncols=4)
        output_dir = cls.params[cf.OUT_DIR]

        for i, t in enumerate(tiles):
            tsincr_file_n = os.path.join(output_dir, 'tsincr_{}.npy'.format(i))
            cls.tsincr_mpi[t.top_left_y:t.bottom_right_y,
                t.top_left_x: t.bottom_right_x, :] = np.load(tsincr_file_n)

            tscum_file_n = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))
            cls.tscum_mpi[t.top_left_y:t.bottom_right_y,
                t.top_left_x: t.bottom_right_x, :] = np.load(tscum_file_n)

        cls.tsincr_tifs_mpi = glob.glob(os.path.join(cls.params[cf.OUT_DIR],
                                                     'tsincr*.tif'))
        cls.tscum_tifs_mpi = glob.glob(os.path.join(cls.params[cf.OUT_DIR],
                                                    'tscum*.tif'))

        cls.tscum_tifs_mpi.sort()
        cls.tsincr_tifs_mpi.sort()

        # remove all temp numpy files after test
        for f in glob.glob(os.path.join(output_dir, '*.npy')):
            os.remove(f)

    def calc_non_mpi_time_series(self):

        # copy sydney_tif files in temp_dir
        self.params[cf.OUT_DIR] = self.temp_dir
        xlks, ylks, crop = run_pyrate.transform_params(self.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            self.base_unw_paths, crop, self.params, xlks)
        # create the dest_paths files
        run_prepifg.gamma_prepifg(self.base_unw_paths, self.params)

        run_pyrate.process_ifgs(dest_paths, self.params)

        # tsvel_file = os.path.join(self.params[cf.OUT_DIR], 'tsvel.npy')
        tsincr_file = os.path.join(self.params[cf.OUT_DIR], 'tsincr.npy')
        tscum_file = os.path.join(self.params[cf.OUT_DIR], 'tscum.npy')
        # self.tsvel = np.load(tsvel_file)
        self.tsincr = np.load(tsincr_file)
        self.tscum = np.load(tscum_file)

        self.tsincr_tifs = glob.glob(os.path.join(self.params[cf.OUT_DIR],
                                                  'tsincr*.tif'))
        self.tscum_tifs = glob.glob(os.path.join(self.params[cf.OUT_DIR],
                                                 'tscum*.tif'))
        self.tsincr_tifs.sort()
        self.tscum_tifs.sort()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir)
        shutil.rmtree(cls.temp_dir)

    def test_mpi_time_series(self):
        for looks, ref_method in product(range(1, 5), [1, 2]):
            print '=======Testing timeseries for looks:', looks, \
                'ref_method: ', ref_method
            self.params[cf.IFG_LKSX] = looks
            self.params[cf.IFG_LKSY] = looks
            self.params[cf.REF_EST_METHOD] = ref_method
            self.process()
            mlooked_ifgs = glob.glob(os.path.join(
                self.tif_dir, '*_{looks}rlks_*cr.tif'.format(looks=looks)))
            self.assertEqual(len(mlooked_ifgs), 17)
            self.calc_non_mpi_time_series()
            np.testing.assert_array_almost_equal(self.tsincr,
                                                 self.tsincr_mpi,
                                                 decimal=4)
            np.testing.assert_array_almost_equal(self.tscum,
                                                 self.tscum_mpi,
                                                 decimal=4)
            # test the tsincr tifs from serial vs mpi from tiles
            for f, g in zip(self.tsincr_tifs, self.tsincr_tifs_mpi):
                fs = gdal.Open(f, gdal.GA_ReadOnly)
                gs = gdal.Open(g, gdal.GA_ReadOnly)
                self.assertDictEqual(fs.GetMetadata(), gs.GetMetadata())
                np.testing.assert_array_almost_equal(fs.ReadAsArray(),
                                                     gs.ReadAsArray(),
                                                     decimal=4)

            # test the tscum from serial vs mpi from tiles
            for f, g in zip(self.tscum_tifs, self.tscum_tifs_mpi):
                fs = gdal.Open(f, gdal.GA_ReadOnly)
                gs = gdal.Open(g, gdal.GA_ReadOnly)
                self.assertDictEqual(fs.GetMetadata(), gs.GetMetadata())
                np.testing.assert_array_almost_equal(fs.ReadAsArray(),
                                                     gs.ReadAsArray(),
                                                     decimal=4)

        log_file = glob.glob(os.path.join(self.tif_dir, '*.log'))[0]
        self.assertTrue(os.path.exists(log_file))


if __name__ == "__main__":
    unittest.main()
