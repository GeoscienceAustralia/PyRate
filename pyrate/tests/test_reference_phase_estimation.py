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
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate.reference_phase_estimation import estimate_ref_phase
from pyrate.scripts import run_prepifg
from pyrate.tests import common
from pyrate import reference_phase_estimation as rpe


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

        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

        ifgs = ifg_instance.ifgs

        for i in ifgs:
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.nodata_value = 0
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        # manually close for windows compatibility
        for i in ifgs:
            i.close()

        cls.matlab_ref_phs = [-18.2191658020020,
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

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(self.matlab_ref_phs, self.ref_phs,
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
        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = True

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

        ifgs = ifg_instance.ifgs

        for i in ifgs:
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.nodata_value = 0
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        # manually close for windows compatibility
        for i in ifgs:
            i.close()

        cls.matlab_ref_phs = [-18.2191658020020,
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

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_out_dir)

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(self.matlab_ref_phs, self.ref_phs,
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

        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

        ifgs = ifg_instance.ifgs

        for c, i in enumerate(ifgs):
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.nodata_value = 0
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        cls.ref_phs_mthod2, cls.ifgs = \
            estimate_ref_phase(ifgs, params, refx, refy)

        cls.matlab_ref_phs = [-21.4459648132324,
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
        np.testing.assert_array_almost_equal(self.matlab_ref_phs,
                                             self.ref_phs_mthod2, decimal=3)


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

        ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

        ifgs = ifg_instance.ifgs

        for c, i in enumerate(ifgs):
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.nodata_value = 0
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        cls.ref_phs_mthod2, cls.ifgs = \
            estimate_ref_phase(ifgs, params, refx, refy)

        cls.matlab_ref_phs = [-21.4459648132324,
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
        np.testing.assert_array_almost_equal(self.matlab_ref_phs,
                                             self.ref_phs_mthod2, decimal=3)


class RefPhaseEstimationMPITest(unittest.TestCase):
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
        cls.params[cf.OUT_DIR] = cls.tif_dir
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.REF_EST_METHOD] = 1
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        cls.base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])

    @classmethod
    def process(cls):
        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            cls.base_unw_paths, crop, cls.params, xlks)

        # create the dest_paths files
        run_prepifg.gamma_prepifg(cls.base_unw_paths, cls.params)

        # now create test ifgs
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

        # give the log file any name
        cls.log_file = os.path.join(cls.tif_dir, 'maxvar_mpi.log')

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

    def calc_non_mpi_time_series(self):
        for i in self.ifgs:
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.nodata_value = 0
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()

        if self.params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(self.ifgs, self.params)

        refx, refy = run_pyrate.find_reference_pixel(self.ifgs, self.params)
        return rpe.estimate_ref_phase(self.ifgs, self.params, refx, refy)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir)

    def test_mpi_mst_single_processor(self):
        # TODO: Why MPI test does not match for looks=2 and ref_phase_method=2
        for looks, ref_method in product([1, 3, 4], [1, 2]):
            self.params[cf.IFG_LKSX] = looks
            self.params[cf.IFG_LKSY] = looks
            self.params[cf.REF_EST_METHOD] = ref_method
            self.process()
            mlooked_ifgs = glob.glob(os.path.join(
                self.tif_dir, '*_{looks}rlks_*cr.tif'.format(looks=looks)))
            self.assertEqual(len(mlooked_ifgs), 17)
            original_ref_phs = self.calc_non_mpi_time_series()[0]
            np.testing.assert_array_almost_equal(original_ref_phs, self.ref_phs,
                                                 decimal=2)

    def test_maxvar_log_written(self):
        self.process()
        log_file = glob.glob(os.path.join(self.tif_dir, '*.log'))[0]
        self.assertTrue(os.path.exists(log_file))


if __name__ == '__main__':
    unittest.main()
