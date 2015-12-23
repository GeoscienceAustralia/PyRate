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
from pyrate.reference_phase_estimation import estimate_ref_phase


class RefPhsEstimationMatlabTest(unittest.TestCase):
    def setUp(self):

        # start each full test run cleanly
        shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

        os.makedirs(SYD_TEST_OUT)

        self.matlab_ref_phs = [-18.2191658020020,
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
        for i in ifgs:
            i.convert_to_mm()
            i.write_modified_phase()

        if params[cf.ORBITAL_FIT] != 0:
            run_pyrate.remove_orbital_error(ifgs, params)

        refx, refy = run_pyrate.find_reference_pixel(ifgs, params)
        ref_phs, ifgs = estimate_ref_phase(ifgs, params, refx, refy)

        np.testing.assert_array_almost_equal(self.matlab_ref_phs, ref_phs,
                                             decimal=4)


if __name__ == '__main__':
    unittest.main()
