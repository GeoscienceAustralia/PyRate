__author__ = 'Sudipta Basak'
__date_created__ = '23/12/15'

import unittest
import os
import sys
import shutil
import numpy as np
import tempfile
import subprocess
import glob
from itertools import product

from pyrate import config as cf
from pyrate.scripts import run_pyrate
from pyrate.tests.common import SYD_TEST_DIR
from pyrate.reference_phase_estimation import (estimate_ref_phase,
                                               ReferencePhaseError)
from pyrate.scripts import run_prepifg
from pyrate.tests import common
from pyrate import ifgconstants as ifc
from pyrate import shared

matlab_ref_phs_method1 = [-18.2191658020020,
                      27.7119445800781,
                      -18.4944229125977,
                      -2.92210483551025,
                      31.1168708801270,
                      21.2123012542725,
                      9.01810073852539,
                      6.08130645751953,
                      -3.79313516616821,
                      -11.3826837539673,
                      -7.28352737426758,
                      17.6365375518799,
                      -12.8673439025879,
                      5.46325922012329,
                      -35.4149475097656,
                      -13.5371961593628,
                      -12.7864856719971]


matlab_ref_phs_method2 = [-21.4459648132324,
                          27.1714553833008,
                          -20.8264484405518,
                          -3.47468209266663,
                          30.4519863128662,
                          22.3201427459717,
                          9.58487224578857,
                          4.81979084014893,
                          -3.89160847663879,
                          -12.0131330490112,
                          -8.64702987670898,
                          19.2060871124268,
                          -9.92049789428711,
                          4.38952684402466,
                          -34.9590339660645,
                          -14.3167810440063,
                          -11.9066228866577]


class RefPhsTests(unittest.TestCase):
    """Basic reference phase estimation tests"""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.params = dict()
        self.params[cf.REF_EST_METHOD] = 1
        self.params[cf.PARALLEL] = False
        self.refpx, self.refpy = 38, 58
        shared.copytree(common.SYD_TEST_TIF, self.tmp_dir)
        sydney_tifs = glob.glob(os.path.join(self.tmp_dir, "*.tif"))
        for s in sydney_tifs:
            os.chmod(s, 0644)
        self.ifgs = common.sydney_data_setup(self.tmp_dir, is_dir=True)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_need_at_least_two_ifgs(self):
        self.assertRaises(ReferencePhaseError, estimate_ref_phase,
                          self.ifgs[:1], self.params, self.refpx, self.refpy)

    def test_metadata(self):
        estimate_ref_phase(self.ifgs, self.params, self.refpx, self.refpy)
        for s in self.ifgs:
            self.assertEqual(s.dataset.GetMetadataItem(ifc.REF_PHASE),
                             ifc.REF_PHASE_REMOVED)

    def test_mixed_metadata_raises(self):
        # correct reference phase for some of the ifgs
        estimate_ref_phase(self.ifgs[:5], self.params, self.refpx, self.refpy)

        # now it should raise exception if we wnat to correct
        # refernece phase again on all of them
        self.assertRaises(ReferencePhaseError, estimate_ref_phase,
                          self.ifgs, self.params, self.refpx, self.refpy)



class RefPhsEstimationMatlabTestMethod1Serial(unittest.TestCase):
    """
    Reference phase estimation method 1 is tested vs matlab output
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

        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = False

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = shared.pre_prepare_ifgs(dest_paths, params)

        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        # end run_pyrate copy

        # manually close for windows compatibility
        for i in ifgs:
            i.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(matlab_ref_phs_method1, self.ref_phs,
                                             decimal=4)

    def test_ifgs_after_reference_phase_estimation(self):
        MATLAB_REF_PHASE_DIR = os.path.join(SYD_TEST_DIR,
                                                     'matlab_ref_phase_est')

        onlyfiles = [f for f in os.listdir(MATLAB_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(MATLAB_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                MATLAB_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data),
                        np.nanmean(self.ifgs[k].phase_data), places=4)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))


class RefPhsEstimationMatlabTestMethod1Parallel(unittest.TestCase):
    """
    Reference phase estimation method 1 is tested vs matlab output
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

        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = True

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = shared.pre_prepare_ifgs(dest_paths, params)

        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        # end run_pyrate copy

        # manually close for windows compatibility
        for i in ifgs:
            i.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(matlab_ref_phs_method1, self.ref_phs,
                                             decimal=4)

    def test_ifgs_after_reference_phase_estimation(self):
        MATLAB_REF_PHASE_DIR = os.path.join(SYD_TEST_DIR,
                                                     'matlab_ref_phase_est')

        onlyfiles = [f for f in os.listdir(MATLAB_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(MATLAB_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                MATLAB_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data),
                        np.nanmean(self.ifgs[k].phase_data), places=4)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))


class RefPhsEstimationMatlabTestMethod2Serial(unittest.TestCase):
    """
    Reference phase estimation method 2 is tested vs matlab output
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

        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = False


        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = shared.pre_prepare_ifgs(dest_paths, params)

        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        # end run_pyrate copy

        for i in ifgs:
            i.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_ifgs_after_reference_phase_estimation(self):
        MATLAB_REF_PHASE_DIR = os.path.join(SYD_TEST_DIR,
                                                     'matlab_ref_phase_est')

        onlyfiles = [f for f in os.listdir(MATLAB_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(MATLAB_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                MATLAB_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data),
                        np.nanmean(self.ifgs[k].phase_data), places=4)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))


        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(matlab_ref_phs_method2,
                                             self.ref_phs, decimal=3)


class RefPhsEstimationMatlabTestMethod2Parallel(unittest.TestCase):
    """
    Reference phase estimation method 2 is tested vs matlab output

    """
    # TODO: Improve the parallel tests to remove duplication from serial tests
    @classmethod
    def setUpClass(cls):

        params = cf.get_config_params(
                os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))

        cls.temp_out_dir = tempfile.mkdtemp()

        sys.argv = ['run_prepifg.py', os.path.join(SYD_TEST_DIR,
                                     'pyrate_system_test.conf')]
        params[cf.OUT_DIR] = cls.temp_out_dir

        run_prepifg.main(params)

        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = True

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = shared.pre_prepare_ifgs(dest_paths, params)
        mst_grid = run_pyrate.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

        # Estimate and remove orbit errors
        run_pyrate.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = shared.pre_prepare_ifgs(dest_paths, params)

        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)
        # end run_pyrate copy

        for i in ifgs:
            i.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_ifgs_after_reference_phase_estimation(self):
        MATLAB_REF_PHASE_DIR = os.path.join(SYD_TEST_DIR,
                                                     'matlab_ref_phase_est')

        onlyfiles = [f for f in os.listdir(MATLAB_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(MATLAB_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                MATLAB_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data),
                        np.nanmean(self.ifgs[k].phase_data), places=4)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))


        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(matlab_ref_phs_method2,
                                             self.ref_phs, decimal=3)


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
        # cls.params[cf.OUT_DIR] = cls.tif_dir
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.APS_CORRECTION] = 0
        cls.params[cf.REF_EST_METHOD] = 1
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        cls.base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])

    @classmethod
    def process(cls):
        cls.params[cf.OUT_DIR] = cls.tif_dir
        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        cls.dest_paths = run_pyrate.get_dest_paths(
            cls.base_unw_paths, crop, cls.params, xlks)

        # create the dest_paths files
        run_prepifg.gamma_prepifg(cls.base_unw_paths, cls.params)

        # now create test ifgs
        cls.ifgs = common.sydney_data_setup(datafiles=cls.dest_paths)
        # give the log file any name
        cls.log_file = os.path.join(cls.tif_dir, 'ref_phs_mpi.log')

        # create the conf file in out_dir
        cls.conf_file = tempfile.mktemp(suffix='.conf', dir=cls.tif_dir)
        cf.write_config_file(cls.params, cls.conf_file)

        # conf file crated?
        assert os.path.exists(cls.conf_file)

        # Calc mst using MPI
        str = 'mpirun -np 2 python pyrate/nci/run_pyrate_pypar.py ' + \
              cls.conf_file
        cmd = str.split()
        subprocess.check_call(cmd)

        # load the ref_phs file for testing
        ref_phs_file = os.path.join(cls.params[cf.OUT_DIR], 'ref_phs.npy')
        cls.ref_phs = np.load(ref_phs_file)

    @classmethod
    def calc_non_mpi_ref_phase(cls):
        cls.tmp_dir = tempfile.mkdtemp()
        cls.params[cf.OUT_DIR] = cls.tmp_dir
        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            cls.base_unw_paths, crop, cls.params, xlks)
        # create the dest_paths files
        run_prepifg.gamma_prepifg(cls.base_unw_paths, cls.params)

        run_pyrate.process_ifgs(dest_paths, cls.params)
        # load the ref_phs file for testing
        ref_phs_file = os.path.join(cls.params[cf.OUT_DIR], 'ref_phs.npy')
        return np.load(ref_phs_file)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir)
        shutil.rmtree(cls.tmp_dir)

    def test_mpi_ref_phase(self):
        for looks, ref_method, orbfit_method in product([1, 2, 3, 4], [1, 2], [1]):
            print 'Testing reference phase looks: {looks}, ' \
                  'ref_method: {ref_method}, ' \
                  'orbfit method: {orbfit}\n'.format(looks=looks,
                                                   ref_method=ref_method,
                                                   orbfit=orbfit_method)
            self.params[cf.IFG_LKSX] = looks
            self.params[cf.IFG_LKSY] = looks
            self.params[cf.REF_EST_METHOD] = ref_method
            self.params[cf.ORBITAL_FIT_METHOD] = orbfit_method
            self.process()
            mlooked_ifgs = glob.glob(os.path.join(
                self.tif_dir, '*_{looks}rlks_*cr.tif'.format(looks=looks)))
            self.assertEqual(len(mlooked_ifgs), 17)
            original_ref_phs = self.calc_non_mpi_ref_phase()
            np.testing.assert_array_almost_equal(original_ref_phs, self.ref_phs,
                                                 decimal=3)

        '''test log generated'''
        log_file = glob.glob(os.path.join(self.tif_dir, '*.log'))[0]
        self.assertTrue(os.path.exists(log_file))


if __name__ == '__main__':
    unittest.main()
