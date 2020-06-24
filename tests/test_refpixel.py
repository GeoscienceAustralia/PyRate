#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains tests for the refpixel.py PyRate module.
"""
import copy
import unittest
import shutil
from subprocess import check_call, run
from pathlib import Path
import pytest
from numpy import nan, mean, std, isnan

from pyrate.core import config as cf
from pyrate.core.refpixel import ref_pixel, _step, RefPixelError
from pyrate.core import shared, ifgconstants as ifc

from pyrate import process
from pyrate.configuration import Configuration
from tests.common import TEST_CONF_ROIPAC
from tests.common import small_data_setup, MockIfg, copy_small_ifg_file_list, \
    copy_and_setup_small_data, manipulate_test_conf, assert_two_dirs_equal, PYTHON3P6


# TODO: figure out how  editing  resource.setrlimit fixes the error
# to fix the open to many files error
# https://stackoverflow.com/questions/18280612/ioerror-errno-24-too-many-open-files

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
        self.ifgs = small_data_setup()
        self.params = cf.get_config_params(TEST_CONF_ROIPAC)
        self.params[cf.REFNX] = REFNX
        self.params[cf.REFNY] = REFNY
        self.params[cf.REF_CHIP_SIZE] = CHIPSIZE
        self.params[cf.REF_MIN_FRAC] = MIN_FRAC
        self.params[cf.PARALLEL] = PARALLEL

    def test_missing_chipsize(self):
        self.params[cf.REF_CHIP_SIZE] = None
        self.assertRaises(cf.ConfigException, ref_pixel, self.ifgs, self.params)

    def test_chipsize_valid(self):
        for illegal in [0, -1, -15, 1, 2, self.ifgs[0].ncols+1, 4, 6, 10, 20]:
            self.params[cf.REF_CHIP_SIZE] = illegal
            self.assertRaises(RefPixelError, ref_pixel, self.ifgs, self.params)

    def test_minimum_fraction_missing(self):
        self.params[cf.REF_MIN_FRAC] = None
        self.assertRaises(cf.ConfigException, ref_pixel, self.ifgs, self.params)

    def test_minimum_fraction_threshold(self):
        for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
            self.params[cf.REF_MIN_FRAC] = illegal
            self.assertRaises(RefPixelError, ref_pixel, self.ifgs, self.params)

    def test_search_windows(self):
        # 45 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 46, 50, 100]:
            self.params[cf.REFNX] = illegal
            self.assertRaises(RefPixelError, ref_pixel, self.ifgs, self.params)

        # 40 is max # cells a width 3 sliding window can iterate over
        for illegal in [-5, -1, 0, 71, 85, 100]:
            self.params[cf.REFNY] = illegal
            self.assertRaises(RefPixelError, ref_pixel, self.ifgs, self.params)

    def test_missing_search_windows(self):
        self.params[cf.REFNX] = None
        self.assertRaises(cf.ConfigException, ref_pixel, self.ifgs, self.params)

        self.params[cf.REFNX] = REFNX
        self.params[cf.REFNY] = None

        self.assertRaises(cf.ConfigException, ref_pixel, self.ifgs, self.params)


class ReferencePixelTests(unittest.TestCase):
    """
    Tests reference pixel search
    """

    def setUp(self):
        self.params = cf.get_config_params(TEST_CONF_ROIPAC)
        self.params[cf.OUT_DIR], self.ifgs = copy_and_setup_small_data()
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
        self.assertRaises(ValueError, ref_pixel, mock_ifgs, self.params)

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


class LegacyEqualityTest(unittest.TestCase):

    def setUp(self):
        self.params = cf.get_config_params(TEST_CONF_ROIPAC)
        self.params[cf.PARALLEL] = 0
        self.params[cf.OUT_DIR], self.ifg_paths = copy_small_ifg_file_list()
        conf_file = Path(self.params[cf.OUT_DIR], 'conf_file.conf')
        cf.write_config_file(params=self.params, output_conf_file=conf_file)
        self.params = Configuration(conf_file).__dict__
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

    def test_small_test_data_ref_pixel_lat_lon_provided(self):
        self.params[cf.REFX], self.params[cf.REFY] = 150.941666654, -34.218333314
        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_small_test_data_ref_pixel(self):
        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_small_test_data_ref_chipsize_15(self):

        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_chipsize_15)
        self.assertEqual(refx, 7)
        self.assertEqual(refy, 7)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_metadata(self):
        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_chipsize_15)
        for i in self.ifg_paths:
            ifg = shared.Ifg(i)
            ifg.open(readonly=True)
            md = ifg.meta_data
            for k, v in zip([ifc.PYRATE_REFPIX_X, ifc.PYRATE_REFPIX_Y, ifc.PYRATE_REFPIX_LAT,
                            ifc.PYRATE_REFPIX_LON, ifc.PYRATE_MEAN_REF_AREA, ifc.PYRATE_STDDEV_REF_AREA],
                            [str(refx), str(refy), 0, 0, 0, 0]):
                assert k in md  # metadata present
                # assert values
            ifg.close()

    def test_small_test_data_ref_all_1(self):
        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_all_1s)
        self.assertAlmostEqual(0.7, self.params_all_1s[cf.REF_MIN_FRAC])
        self.assertEqual(1, self.params_all_1s[cf.REFNX])
        self.assertEqual(1, self.params_all_1s[cf.REFNY])
        self.assertEqual(refx, 2)
        self.assertEqual(refy, 2)


class LegacyEqualityTestMultiprocessParallel(unittest.TestCase):

    def setUp(self):
        self.params = cf.get_config_params(TEST_CONF_ROIPAC)
        self.params[cf.PARALLEL] = True
        self.params[cf.OUT_DIR], self.ifg_paths = copy_small_ifg_file_list()
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

    def test_small_test_data_ref_pixel(self):
        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.8, self.params[cf.REF_MIN_FRAC])

    def test_more_small_test_data_ref_pixel(self):

        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_alt_ref_frac)
        self.assertEqual(refx, 38)
        self.assertEqual(refy, 58)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_small_test_data_ref_pixel_all_2(self):

        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_all_2s)
        self.assertEqual(refx, 25)
        self.assertEqual(refy, 2)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_small_test_data_ref_chipsize_15(self):

        refx, refy = process._ref_pixel_calc(self.ifg_paths,
                                             self.params_chipsize_15)
        self.assertEqual(refx, 7)
        self.assertEqual(refy, 7)
        self.assertAlmostEqual(0.5, self.params_alt_ref_frac[cf.REF_MIN_FRAC])

    def test_small_test_data_ref_all_1(self):

        refx, refy = process._ref_pixel_calc(self.ifg_paths, self.params_all_1s)

        self.assertAlmostEqual(0.7, self.params_all_1s[cf.REF_MIN_FRAC])
        self.assertEqual(1, self.params_all_1s[cf.REFNX])
        self.assertEqual(1, self.params_all_1s[cf.REFNY])
        self.assertEqual(refx, 2)
        self.assertEqual(refy, 2)


@pytest.mark.slow
@pytest.mark.skip(PYTHON3P6, reason='Skipped in python3p6')
def test_error_msg_refpixel_out_out_bounds(tempdir, gamma_conf):
    "check correct latitude/longitude refpixel error is raised when specified refpixel is out of bounds"
    for x, (refx, refy) in zip(['longitude', 'latitude', 'longitude and latitude'],
                               [(150., -34.218333314), (150.941666654, -34.), (150, -34)]):
        _, out = _get_mlooked_files(gamma_conf, Path(tempdir()), refx=refx, refy=refy)
        msg = "Supplied {} value is outside the bounds of the interferogram data"
        assert msg.format(x) in out.stderr


@pytest.mark.slow
@pytest.mark.skip(PYTHON3P6, reason='Skipped in python3p6')
def test_gamma_ref_pixel_search_vs_lat_lon(tempdir, gamma_conf):
    params_1, _ = _get_mlooked_files(gamma_conf, Path(tempdir()), refx=-1, refy=-1)
    params_2, _ = _get_mlooked_files(gamma_conf, Path(tempdir()), refx=150.941666654, refy=-34.218333314)
    assert_two_dirs_equal(params_1[cf.OUT_DIR], params_2[cf.OUT_DIR], f"*{params_1[cf.IFG_CROP_OPT]}cr.tif", 18)


def _get_mlooked_files(gamma_conf, tdir, refx, refy):
    params = manipulate_test_conf(gamma_conf, tdir)
    params[cf.REFX] = refx
    params[cf.REFY] = refy
    output_conf_file = 'config.conf'
    output_conf = tdir.joinpath(output_conf_file)
    cf.write_config_file(params=params, output_conf_file=output_conf)
    check_call(f"pyrate conv2tif -f {output_conf}", shell=True)
    check_call(f"pyrate prepifg -f {output_conf}", shell=True)
    stdout = run(f"pyrate process -f {output_conf}", shell=True, capture_output=True, text=True)
    print("============================================", stdout)
    return params, stdout
