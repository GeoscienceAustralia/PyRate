__author__ = 'Sudipta Basak'
__date_created__ = '23/12/15'

import unittest
import os
import sys
import shutil
import numpy as np
import tempfile

from pyrate import config as cf
from pyrate.scripts import run_pyrate
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate.reference_phase_estimation import estimate_ref_phase
from pyrate.scripts import run_prepifg
from pyrate.tests import common


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


class RefPhsEstimationMatlabTestMethod2(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()
