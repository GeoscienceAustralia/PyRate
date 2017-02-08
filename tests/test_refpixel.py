"""
Collection of tests for validating PyRate's reference pixel code.
"""
import copy
import os
import unittest
import pytest
import tempfile
import shutil
from numpy import nan, mean, std, isnan

from pyrate import config as cf
from pyrate.config import ConfigException
from pyrate.refpixel import ref_pixel, RefPixelError, step
from pyrate.scripts import run_prepifg
from pyrate.scripts import run_pyrate
from pyrate import mpiops
from tests import common
from tests.common import SYD_TEST_DIR
from tests.common import sydney_data_setup, MockIfg, sydney_ifg_file_list

# default testing values
REFNX = 5
REFNY = 7
MIN_FRAC = 0.7
CHIPSIZE = 3
PARALLEL = False


class ReferencePixelInputTests(unittest.TestCase):
    '''
    Verifies error checking capabilities of the reference pixel function
    '''

    def setUp(self):
        self.ifgs = sydney_data_setup()
        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))
        self.params[cf.REFNX] = REFNX
        self.params[cf.REFNY] = REFNY
        self.params[cf.REF_CHIP_SIZE] = CHIPSIZE
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL

    def test_missing_chipsize(self):
        self.params[cf.REF_CHIP_SIZE] = None
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, self.params)

    def test_chipsize_valid(self):
        for illegal in [0, -1, -15, 1, 2, self.ifgs[0].ncols+1, 4, 6, 10, 20]:
            self.params[cf.REF_CHIP_SIZE] = illegal
            self.assertRaises(ValueError, ref_pixel, self.ifgs, self.params)

    def test_minimum_fraction_missing(self):
        self.params[cf.REF_MIN_FRAC] = None
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, self.params)

    def test_minimum_fraction_threshold(self):
        for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
            self.params[cf.REF_MIN_FRAC] = illegal
            self.assertRaises(ValueError, ref_pixel, self.ifgs, self.params)

    def test_search_windows(self):
        # 45 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 46, 50, 100]:
            self.params[cf.REFNX] = illegal
            self.assertRaises(ValueError, ref_pixel, self.ifgs, self.params)

        # 40 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 71, 85, 100]:
            self.params[cf.REFNY] = illegal
            self.assertRaises(ValueError, ref_pixel, self.ifgs, self.params)

    def test_missing_search_windows(self):
        self.params[cf.REFNX] = None
        self.assertRaises(ConfigException, ref_pixel, self.ifgs, self.params)

        self.params[cf.REFNX] = REFNX
        self.params[cf.REFNY] = None

        self.assertRaises(ConfigException, ref_pixel, self.ifgs, self.params)


class ReferencePixelTests(unittest.TestCase):
    """
    Tests reference pixel search
    """

    def setUp(self):
        self.ifgs = sydney_data_setup()
        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))
        self.params[cf.REFNX] = REFNX
        self.params[cf.REFNY] = REFNY
        self.params[cf.REF_CHIP_SIZE] = CHIPSIZE
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL

    def test_all_below_threshold_exception(self):
        # test failure when no valid stacks in dataset

        # rig mock data to be below threshold
        mock_ifgs = [MockIfg(i, 6, 7) for i in self.ifgs]
        for m in mock_ifgs:
            m.phase_data[:1] = nan
            m.phase_data[1:5] = 0.1
            m.phase_data[5:] = nan

        self.params[cf.REFNX] = 2
        self.params[cf.REFNY] = 2
        self.params[cf.REF_CHIP_SIZE] = CHIPSIZE
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL
        self.assertRaises(RefPixelError, ref_pixel, mock_ifgs, self.params)

    def test_refnxy_step_1(self):
        # test step of 1 for refnx|y gets the reference pixel for axis centre
        mock_ifgs = [MockIfg(i, 47, 72) for i in self.ifgs]
        for m in mock_ifgs:
            m.phase_data[:1] = 0.2
            m.phase_data[1:5] = 0.1
            m.phase_data[5:] = 0.3
        exp_refpx = (1, 1)
        self.params[cf.REFNX] = 1
        self.params[cf.REFNY] = 1
        self.params[cf.REF_CHIP_SIZE] = CHIPSIZE
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL
        res = ref_pixel(mock_ifgs, self.params)
        self.assertEqual(exp_refpx, res)

    def test_large_window(self):
        # 5x5 view over a 5x5 ifg with 1 window/ref pix search
        chps = 5
        mockifgs = [MockIfg(i, chps, chps) for i in self.ifgs]
        self.params[cf.REFNX] = 1
        self.params[cf.REFNY] = 1
        self.params[cf.REF_CHIP_SIZE] = chps
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL
        res = ref_pixel(mockifgs, self.params)
        self.assertEqual((2, 2), res)

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
        act = step(width, refnx, radius)
        assert_equal(act, exp)

        # test with 3 windows
        refnx = 3
        exp = [2, 17, 32]
        act = step(width, refnx, radius)
        assert_equal(act, exp)

        # test 4 search windows
        refnx = 4
        exp = [2, 13, 24, 35]
        act = step(width, refnx, radius)
        assert_equal(act, exp)

    def test_ref_pixel(self):
        exp_refpx = (2, 25)
        self.params[cf.REFNX] = 2
        self.params[cf.REFNY] = 2
        self.params[cf.REF_CHIP_SIZE] = 5
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL
        res = ref_pixel(self.ifgs, self.params)
        self.assertEqual(res, exp_refpx)

        # Invalidate first data stack, get new refpix coods & retest
        for i in self.ifgs:
            i.phase_data[:30, :50] = nan

        exp_refpx = (38, 2)
        res = ref_pixel(self.ifgs, self.params)
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
    assert isnan([ulm, urm, llm, lrm]).any() is False

    # coords of the smallest mean is the result
    mn = [ulm, urm, llm, lrm]


class MatlabEqualityTest(unittest.TestCase):

    def setUp(self):
        self.ifg_paths = sydney_ifg_file_list()
        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))
        self.params[cf.PARALLEL] = False
        self.params[cf.OUT_DIR] = tempfile.mkdtemp()
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

    def tearDown(self):
        shutil.rmtree(self.params[cf.OUT_DIR])

    def test_sydney_test_data_ref_pixel(self):
        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_more_sydney_test_data_ref_pixel(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_alt_ref_frac)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_pixel_all_2(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_all_2s)
        self.assertEqual(refx, 25)
        self.assertEqual(refy, 2)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_chipsize_15(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_chipsize_15)
        self.assertEqual(refx, 7)
        self.assertEqual(refy, 7)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_all_1(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_all_1s)

        self.assertAlmostEqual(0.7, self.params_all_1s[cf.REF_MIN_FRAC])
        self.assertEqual(1, self.params_all_1s[cf.REFNX])
        self.assertEqual(1, self.params_all_1s[cf.REFNY])
        self.assertEqual(refx, 2)
        self.assertEqual(refy, 2)


class MatlabEqualityTestMultiprocessParallel(unittest.TestCase):

    def setUp(self):
        self.ifg_paths = sydney_ifg_file_list()
        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))
        self.params[cf.PARALLEL] = True
        self.params[cf.OUT_DIR] = tempfile.mkdtemp()

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

    def tearDown(self):
        shutil.rmtree(self.params[cf.OUT_DIR])

    def test_sydney_test_data_ref_pixel(self):
        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_more_sydney_test_data_ref_pixel(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_alt_ref_frac)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_pixel_all_2(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_all_2s)
        self.assertEqual(refx, 25)
        self.assertEqual(refy, 2)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_chipsize_15(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_chipsize_15)
        self.assertEqual(refx, 7)
        self.assertEqual(refy, 7)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_sydney_test_data_ref_all_1(self):

        refx, refy = run_pyrate.ref_pixel_calc(self.ifg_paths,
                                               self.params_all_1s)

        self.assertAlmostEqual(0.7, self.params_all_1s[cf.REF_MIN_FRAC])
        self.assertEqual(1, self.params_all_1s[cf.REFNX])
        self.assertEqual(1, self.params_all_1s[cf.REFNY])
        self.assertEqual(refx, 2)
        self.assertEqual(refy, 2)


@pytest.fixture(params=range(1, 6))
def modify_config(request, get_config):
    test_conf = common.SYDNEY_TEST_CONF
    params_dict = get_config(test_conf)
    params_dict[cf.IFG_LKSX] = request.param
    params_dict[cf.IFG_LKSY] = request.param
    params_dict[cf.IFG_FILE_LIST] = os.path.join(
        common.SYD_TEST_GAMMA, 'ifms_17')
    params_dict[cf.PARALLEL] = 0
    params_dict[cf.APS_CORRECTION] = 0
    return params_dict


def test_refpixel_mpi(mpisync, tempdir, modify_config,
                      ref_est_method, roipac_or_gamma):
    params_dict = modify_config

    # use the same output dir for all processes
    outdir = mpiops.run_once(tempdir)
    params_dict[cf.OUT_DIR] = outdir
    params_dict[cf.REF_EST_METHOD] = ref_est_method
    xlks, ylks, crop = cf.transform_params(params_dict)
    if roipac_or_gamma == 0:
        params_dict[cf.IFG_FILE_LIST] = os.path.join(
            common.SYD_TEST_OBS, 'ifms_17')

    base_unw_paths = cf.original_ifg_paths(params_dict[cf.IFG_FILE_LIST])
    # dest_paths are tifs that have been geotif converted and multilooked

    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params_dict, xlks)
    # run prepifg, create the dest_paths files
    if mpiops.rank == 0:
        # create the dest_paths files
        if roipac_or_gamma == 0:
            run_prepifg.roipac_prepifg(base_unw_paths, params_dict)
        else:
            run_prepifg.gamma_prepifg(base_unw_paths, params_dict)

    mpiops.comm.barrier()

    refpx, refpy = run_pyrate.ref_pixel_calc(dest_paths, params_dict)

    # old single process ref pixel calc
    if mpiops.rank == 0:
        ifgs = sydney_data_setup(datafiles=dest_paths)
        refy, refx = ref_pixel(ifgs, params_dict)
        for i in ifgs:  # close ifgs
            i.close()
        assert (refpx, refpy) == (refx, refy)
        shutil.rmtree(outdir)


if __name__ == "__main__":
    unittest.main()
