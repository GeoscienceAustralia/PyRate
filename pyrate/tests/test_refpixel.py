'''
Collection of tests for validating PyRate's reference pixel code.

.. codeauthor:: Ben Davies
'''

import unittest
import os
import copy
from numpy import nan, mean, std, isnan

from common import sydney_data_setup, MockIfg
from pyrate.refpixel import ref_pixel, RefPixelError, _step
from pyrate.config import ConfigException
from pyrate.scripts import run_pyrate
from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR
from pyrate import config as cf

# default testing values
REFNX = 5
REFNY = 7
MIN_FRAC = 0.7
CHIPSIZE = 3


class ReferencePixelInputTests(unittest.TestCase):
    '''
    Verifies error checking capabilities of the reference pixel function
    '''

    def setUp(self):
        self.ifgs = sydney_data_setup()

    def test_missing_chipsize(self):
        args = (REFNX, REFNY, None, MIN_FRAC)
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, *args)


    def test_chipsize_valid(self):
        for illegal in [0, -1, -15, 1, 2, self.ifgs[0].ncols+1, 4, 6, 10, 20]:
            args = (REFNX, REFNY, illegal, MIN_FRAC)
            self.assertRaises(ValueError, ref_pixel, self.ifgs, *args)


    def test_minimum_fraction_missing(self):
        args = (REFNX, REFNY, CHIPSIZE, None)
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, *args)


    def test_minimum_fraction_threshold(self):
        for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
            args = (REFNX, REFNY, CHIPSIZE, illegal)
            self.assertRaises(ValueError, ref_pixel, self.ifgs, *args)


    def test_search_windows(self):
        # 45 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 46, 50, 100]:
            args = (illegal, REFNY, CHIPSIZE, MIN_FRAC)
            self.assertRaises(ValueError, ref_pixel, self.ifgs, *args)

        # 40 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 71, 85, 100]:
            args = (REFNX, illegal, CHIPSIZE, MIN_FRAC)
            self.assertRaises(ValueError, ref_pixel, self.ifgs, *args)


    def test_missing_search_windows(self):
        args = (None, REFNY, CHIPSIZE, MIN_FRAC)
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, *args)

        args = (REFNX, None, CHIPSIZE, MIN_FRAC)
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, *args)



class ReferencePixelTests(unittest.TestCase):
    """
    Tests reference pixel search
    """

    def setUp(self):
        self.ifgs = sydney_data_setup()

    def test_all_below_threshold_exception(self):
        # test failure when no valid stacks in dataset

        # rig mock data to be below threshold
        mock_ifgs = [MockIfg(i, 6, 7) for i in self.ifgs]
        for m in mock_ifgs:
            m.phase_data[:1] = nan
            m.phase_data[1:5] = 0.1
            m.phase_data[5:] = nan

        args = (2, 2, 3, 0.7)
        self.assertRaises(RefPixelError, ref_pixel, mock_ifgs, *args)

    def test_refnxy_step_1(self):
        # test step of 1 for refnx|y gets the reference pixel for axis centre
        mock_ifgs = [MockIfg(i, 47, 72) for i in self.ifgs]
        for m in mock_ifgs:
            m.phase_data[:1] = 0.2
            m.phase_data[1:5] = 0.1
            m.phase_data[5:] = 0.3
        exp_refpx = (1, 1)
        res = ref_pixel(mock_ifgs, refnx=1, refny=1, chipsize=3, min_frac=0.7)
        self.assertEqual(exp_refpx, res)

    def test_large_window(self):
        # 5x5 view over a 5x5 ifg with 1 window/ref pix search
        chps = 5
        mockifgs = [MockIfg(i, chps, chps) for i in self.ifgs]
        res = ref_pixel(mockifgs, refnx=1, refny=1, chipsize=chps, min_frac=0.7)
        self.assertEqual((2,2), res)

    def test_step(self):
        # test different search windows to verify x/y step calculation

        # convenience testing function
        def assert_equal(actual, expected):
            for a, e in zip(actual, expected):
                self.assertEqual(a, e)

        # start with simple corner only test
        width = 47
        radius = 2
        refnx = 2
        exp = [2, 25, 44]
        act = _step(width, refnx, radius)
        assert_equal(act, exp)

        # test with 3 windows
        refnx = 3
        exp = [2, 17, 32]
        act = _step(width, refnx, radius)
        assert_equal(act, exp)

        # test 4 search windows
        refnx = 4
        exp = [2, 13, 24, 35]
        act = _step(width, refnx, radius)
        assert_equal(act, exp)

    def test_ref_pixel(self):
        exp_refpx = (2, 25)
        res = ref_pixel(self.ifgs, 2, 2, 5, 0.7)
        self.assertEqual(res, exp_refpx)

        # Invalidate first data stack, get new refpix coods & retest
        for i in self.ifgs:
            i.phase_data[:30, :50] = nan

        exp_refpx = (38, 2)
        res = ref_pixel(self.ifgs, 2, 2, 5, 0.7)
        self.assertEqual(res, exp_refpx)


def _expected_ref_pixel(ifgs, cs):
    """Helper function for finding reference pixel when refnx/y=2"""

    # calculate expected data
    data = [i.phase_data for i in ifgs]  # len 17 list of arrays
    ul = [i[:cs, :cs] for i in data]  # upper left corner stack
    ur = [i[:cs, -cs:] for i in data]
    ll = [i[-cs:, :cs] for i in data]
    lr = [i[-cs:, -cs:] for i in data]

    ulm = mean([std(i[~isnan(i)]) for i in ul])  # mean std of all the layers
    urm = mean([std(i[~isnan(i)]) for i in ur])
    llm = mean([std(i[~isnan(i)]) for i in ll])
    lrm = mean([std(i[~isnan(i)]) for i in lr])
    assert isnan([ulm, urm, llm, lrm]).any() == False

    # coords of the smallest mean is the result
    mn = [ulm, urm, llm, lrm]
    print mn, min(mn), mn.index(min(mn))


class MatlabEqualityTest(unittest.TestCase):

    def setUp(self):
        self.ifgs = sydney_data_setup()
        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))

        self.params_alt_ref_frac = copy.copy(self.params)
        self.params_alt_ref_frac[cf.REF_MIN_FRAC] = 0.5
        self.params_all_2s = copy.copy(self.params)
        self.params_all_2s[cf.REFNX] = 2
        self.params_all_2s[cf.REFNY] = 2
        self.params_chipsize_15 = copy.copy(self.params_all_2s)
        self.params_chipsize_15[cf.REF_CHIP_SIZE] = 15
        self.params_all_1s = copy.copy(self.params)
        self.params_all_1s[cf.REFNX] = 1
        self.params_all_1s[cf.REFNY] = 1
        self.params_all_1s[cf.REF_MIN_FRAC] = 0.7

    def test_sydney_test_data_ref_pixel(self):
        refx, refy = run_pyrate.find_reference_pixel(self.ifgs, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_more_sydney_test_data_ref_pixel(self):

        refx, refy = run_pyrate.find_reference_pixel(self.ifgs,
                                                     self.params_alt_ref_frac)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_pixel_all_2(self):

        refx, refy = run_pyrate.find_reference_pixel(self.ifgs,
                                                     self.params_all_2s)
        self.assertEqual(refx, 25)
        self.assertEqual(refy, 2)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_chipsize_15(self):

        refx, refy = run_pyrate.find_reference_pixel(self.ifgs,
                                                     self.params_chipsize_15)
        self.assertEqual(refx, 7)
        self.assertEqual(refy, 7)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_all_1(self):

        refx, refy = run_pyrate.find_reference_pixel(self.ifgs,
                                                     self.params_all_1s)

        self.assertAlmostEqual(0.7, self.params_all_1s[cf.REF_MIN_FRAC])
        self.assertEqual(1, self.params_all_1s[cf.REFNX])
        self.assertEqual(1, self.params_all_1s[cf.REFNY])
        self.assertEqual(refx, 2)
        self.assertEqual(refy, 2)


if __name__ == "__main__":
    unittest.main()
