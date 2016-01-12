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

from pyrate.scripts import run_pyrate, run_prepifg
from pyrate import matlab_mst_kruskal as matlab_mst
from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT
from pyrate.tests.common import SYD_TEST_DIR
from pyrate import config as cf
from pyrate import reference_phase_estimation as rpe
from pyrate import vcm
from pyrate.config import LR_PTHRESH, LR_MAXSIG, LR_NSIG
from pyrate.linrate import linear_rate
from pyrate import mst
from pyrate.vcm import cvd, get_vcmt
from pyrate.tests.common import sydney_data_setup, sydney_data_setup_ifg_file_list

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

#def default_params():
#    return { LR_PTHRESH : 10, LR_NSIG : 3, LR_MAXSIG : 2 }
LR_PTHRESH = 3
LR_NSIG = 3
LR_MAXSIG = 2


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
        vcm = eye(6, 6)
        mst = ones((6, 1, 1))
        mst[4] = 0
        rate, error, samples = linear_rate(self.ifgs, vcm, LR_PTHRESH,
                                           LR_NSIG, LR_MAXSIG, mst)
        assert_array_almost_equal(rate, exprate)
        assert_array_almost_equal(error, experr)
        assert_array_almost_equal(samples, expsamp)
        

class MatlabEqualityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # start each full test run cleanly
        shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

        os.makedirs(SYD_TEST_OUT)

        params = cf.get_config_params(
                os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))
        params[cf.REF_EST_METHOD] = 2

        sys.argv = ['run_prepifg.py', os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR,
                                     'orbital_error.conf')]
        run_prepifg.main()

        xlks, ylks, crop = run_pyrate.transform_params(params)

        base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params, xlks)

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

        # Calculate linear rate map
        params[cf.PARALLEL] = 1
        cls.rate, cls.error, cls.samples = run_pyrate.calculate_linear_rate(
            ifgs, params, vcmt, mst=mst_grid)

        params[cf.PARALLEL] = 2
        cls.rate_2, cls.error_2, cls.samples_2 = \
            run_pyrate.calculate_linear_rate(ifgs, params, vcmt, mst=mst_grid)

        MATLAB_LINRATE_DIR = os.path.join(SYD_TEST_DIR, 'matlab_linrate')

        cls.rate_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'stackmap.csv'), delimiter=',')
        cls.error_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'errormap.csv'), delimiter=',')

        cls.samples_matlab = np.genfromtxt(
            os.path.join(MATLAB_LINRATE_DIR, 'coh_sta.csv'), delimiter=',')

        params[cf.PARALLEL] = 0
        # Calculate linear rate map
        cls.rate_s, cls.error_s, cls.samples_s = \
            run_pyrate.calculate_linear_rate(ifgs, params, vcmt, mst=mst_grid)

    def test_linear_rate_full_parallel(self):
        np.testing.assert_array_almost_equal(
            self.rate[:11, :45], self.rate_matlab[:11, :45], decimal=4)

    def test_lin_rate_error_parallel(self):
        np.testing.assert_array_almost_equal(
            self.error[:11, :45], self.error_matlab[:11, :45], decimal=4)

    def test_lin_rate_samples_parallel(self):
        np.testing.assert_array_almost_equal(
            self.samples[:11, :45], self.samples_matlab[:11, :45], decimal=4)

    def test_linear_rate_full_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.rate_2[:11, :45], self.rate_matlab[:11, :45], decimal=4)

    def test_lin_rate_error_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.error_2[:11, :45], self.error_matlab[:11, :45], decimal=4)

    def test_lin_rate_samples_parallel_pixel_level(self):
        np.testing.assert_array_almost_equal(
            self.samples_2[:11, :45], self.samples_matlab[:11, :45], decimal=4)

    def test_linear_rate_full_serial(self):
        np.testing.assert_array_almost_equal(
            self.rate_s[:11, :45], self.rate_matlab[:11, :45], decimal=4)

    def test_lin_rate_error_serial(self):
        np.testing.assert_array_almost_equal(
            self.error_s[:11, :45], self.error_matlab[:11, :45], decimal=4)

    def test_lin_rate_samples_serial(self):
        np.testing.assert_array_almost_equal(
            self.samples_s[:11, :45], self.samples_matlab[:11, :45], decimal=4)

    # TODO: investigae why the whole arrays don't equal

if __name__ == "__main__":
    unittest.main()
