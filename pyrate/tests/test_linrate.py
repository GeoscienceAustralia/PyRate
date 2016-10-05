'''
CUnittests for linrate.py

.. codeauthor:: Ben Davies, Sudipta Basak
'''

import unittest

from numpy import eye, array, ones, where
from numpy.testing import assert_array_almost_equal
import os
import shutil
import sys
import numpy as np
import tempfile
import subprocess
import glob
from itertools import product

from pyrate.scripts import run_pyrate, run_prepifg
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate import config as cf
from pyrate import reference_phase_estimation as rpe
from pyrate import vcm
from pyrate.config import LR_PTHRESH, LR_MAXSIG, LR_NSIG
from pyrate.linrate import linear_rate
from pyrate.tests import common
from pyrate import vcm as vcm_module
from pyrate import prepifg
from pyrate import shared

# TODO: linear rate code
# 1. replace MST key:value date:date pairs with lists of Ifgs?
# 2.: MG: fake some data for a 1 pixel test of stack()
#     figure out what partial inputs we need to do this as a simple test
# 3. MG is going to check if someome at GA has coded lscov in Py (or C/C++/Fortran??)
#
# MG thinks we'll need:
# 1 all the func args
# 2 a VCM
# 3 the time spans (already in Ifgs?)
# 4 the ifg displacement obs (phase_data)



# class LinearRateTests(unittest.TestCase):
#
#     def test_stack_basic(self):
#         raise NotImplementedError
#
#
#     def test_args(self):
#         raise NotImplementedError("Need sanity tests for args to stack()")

def default_params():
    return {'pthr': 3, 'nsig': 3, 'maxsig': 2, 'parallel': 1, 'processes': 8}


class SinglePixelIfg(object):

    def __init__(self,timespan,phase):
        self.time_span = timespan
        self.phase_data = array([[phase]])


class LinearRateTests(unittest.TestCase):
    """Tests the weighted least squares algorithm for determinining
    the best fitting velocity"""

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

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params, xlks)

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

        # Calculate linear rate map
        params[cf.PARALLEL] = 1
        cls.rate, cls.error, cls.samples = run_pyrate.calculate_linear_rate(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        cls.rate_2, cls.error_2, cls.samples_2 = \
            run_pyrate.calculate_linear_rate(ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 0
        # Calculate linear rate map
        cls.rate_s, cls.error_s, cls.samples_s = \
            run_pyrate.calculate_linear_rate(ifgs, params, vcmt, mst=mst_grid)

        MATLAB_LINRATE_DIR = os.path.join(SYD_TEST_DIR, 'matlab_linrate')

        cls.rate_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'stackmap.csv'), delimiter=',')
        cls.error_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'errormap.csv'), delimiter=',')

        cls.samples_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'coh_sta.csv'), delimiter=',')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_linear_rate_full_parallel(self):
        np.testing.assert_array_almost_equal(
            self.rate, self.rate_matlab, decimal=3)

    def test_lin_rate_error_parallel(self):
        np.testing.assert_array_almost_equal(
            self.error, self.error_matlab, decimal=3)

    def test_lin_rate_samples_parallel(self):
        np.testing.assert_array_almost_equal(
            self.samples, self.samples_matlab, decimal=3)

    def test_linear_rate_full_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.rate_2, self.rate_matlab, decimal=3)

    def test_lin_rate_error_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.error_2, self.error_matlab, decimal=3)

    def test_lin_rate_samples_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.samples_2, self.samples_matlab, decimal=3)

    def test_linear_rate_full_serial(self):
        np.testing.assert_array_almost_equal(
            self.rate_s, self.rate_matlab, decimal=3)

    def test_lin_rate_error_serial(self):
        np.testing.assert_array_almost_equal(
            self.error_s, self.error_matlab, decimal=3)

    def test_lin_rate_samples_serial(self):
        np.testing.assert_array_almost_equal(
            self.samples_s, self.samples_matlab, decimal=3)


class MPITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        cls.params[cf.IFG_FILE_LIST] = os.path.join(
            common.SYD_TEST_GAMMA, 'ifms_17')
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.REF_EST_METHOD] = 1
        cls.params[cf.APS_CORRECTION] = 0
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
        cls.log_file = os.path.join(cls.tif_dir, 'linrate_mpi.log')

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

        rate_file = os.path.join(cls.params[cf.OUT_DIR], 'linrate.npy')
        cls.rate_mpi = np.load(rate_file)
        error_file = os.path.join(cls.params[cf.OUT_DIR], 'linerror.npy')
        cls.error_mpi = np.load(error_file)
        samples_file = os.path.join(cls.params[cf.OUT_DIR], 'linsamples.npy')
        cls.samples_mpi = np.load(samples_file)

    def calc_non_mpi_time_series(self):
        temp_dir = tempfile.mkdtemp()
        # copy sydney_tif files in temp_dir
        self.params[cf.OUT_DIR] = temp_dir
        xlks, ylks, crop = run_pyrate.transform_params(self.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            self.base_unw_paths, crop, self.params, xlks)
        # create the dest_paths files
        run_prepifg.gamma_prepifg(self.base_unw_paths, self.params)

        _, _, _, _, self.rate, self.error, self.samples = \
            run_pyrate.process_ifgs(dest_paths, self.params)

        shutil.rmtree(temp_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir)

    def test_mpi_mst_single_processor(self):
        for looks, ref_method, orbfitmethod in product(range(1, 5), [1, 2],
                                                       [1]):
            print '=======Testing linear rate for looks:', looks, \
                'ref_method: ', ref_method
            self.params[cf.IFG_LKSX] = looks
            self.params[cf.IFG_LKSY] = looks
            self.params[cf.REF_EST_METHOD] = ref_method
            self.params[cf.ORBITAL_FIT_METHOD] = orbfitmethod
            self.process()
            mlooked_ifgs = glob.glob(os.path.join(
                self.tif_dir, '*_{looks}rlks_*cr.tif'.format(looks=looks)))
            self.assertEqual(len(mlooked_ifgs), 17)
            self.calc_non_mpi_time_series()
            np.testing.assert_array_almost_equal(self.rate,
                                                 self.rate_mpi,
                                                 decimal=4)

            np.testing.assert_array_almost_equal(self.error,
                                                 self.error_mpi,
                                                 decimal=4)
            np.testing.assert_array_almost_equal(self.samples,
                                                 self.samples_mpi,
                                                 decimal=4)

        # test logfile was created
        log_file = glob.glob(os.path.join(self.tif_dir, '*.log'))[0]
        self.assertTrue(os.path.exists(log_file))

        # test linerror.tif and linrate.tif was created
        rate_file = glob.glob(os.path.join(self.tif_dir, 'linrate.tif'))[0]
        self.assertTrue(os.path.exists(rate_file))

        error_file = glob.glob(os.path.join(self.tif_dir, 'linerror.tif'))[0]
        self.assertTrue(os.path.exists(error_file))



if __name__ == "__main__":
    unittest.main()
