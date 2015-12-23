__author__ = 'Sudipta Basak'
__date_created__ = '23/12/15'

import unittest
import os
import shutil
from subprocess import call
import numpy as np

from pyrate import config as cf
from pyrate.scripts import run_pyrate
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate.reference_phase_estimation import estimate_ref_phase


class RefPhsEstimationMatlabTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # start each full test run cleanly
        shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

        os.makedirs(SYD_TEST_OUT)
        params = cf.get_config_params(
                os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))

        call(["python", "pyrate/scripts/run_prepifg.py",
              os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf')])

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
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        cls.ref_phs, cls.ifgs = estimate_ref_phase(ifgs, params, refx, refy)

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

if __name__ == '__main__':
    unittest.main()
