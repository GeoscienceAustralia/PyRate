#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
This Python module contains tests for the prepifg.py PyRate module.
"""
import glob
import os
import shutil
import sys
import tempfile
import unittest
from math import floor
from os.path import exists, join

import numpy as np
from numpy import isnan, nanmax, nanmin, ones, nan, reshape, sum as npsum
from numpy.testing import assert_array_almost_equal, assert_array_equal

try:
    from scipy.stats.stats import nanmean
except:  # fix for scipy v0.18.0
    from scipy import nanmean

from osgeo import gdal

from pyrate.scripts import run_prepifg
from pyrate import config as cf
from pyrate.config import mlooked_path
from pyrate.shared import Ifg, DEM
from pyrate.prepifg import CUSTOM_CROP, MAXIMUM_CROP, MINIMUM_CROP, \
    ALREADY_SAME_SIZE
from pyrate.prepifg import prepare_ifgs, resample, PreprocessError, CustomExts
from pyrate.prepifg import extents_from_params
from pyrate.tasks.utils import DUMMY_SECTION_NAME
from pyrate.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR,
    LUIGI,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    NO_DATA_AVERAGING_THRESHOLD,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
    APS_METHOD,
    APS_CORRECTION)

from tests.common import SYD_TEST_MATLAB_PREPIFG_DIR
from tests.common import PREP_TEST_TIF, SYD_TEST_DEM_DIR
from tests.common import SYD_TEST_DEM_TIF
from tests import common

gdal.UseExceptions()

if not exists(PREP_TEST_TIF):
    sys.exit("ERROR: Missing 'prepifg' dir for unittests\n")


# convenience ifg creation funcs
def diff_exts_ifgs():
    """Returns pair of test Ifgs with different extents"""
    bases = ['geo_060619-061002.tif', 'geo_070326-070917.tif']
    random_dir = tempfile.mkdtemp()
    for p in bases:
        shutil.copy(src=os.path.join(PREP_TEST_TIF, p),
                    dst=os.path.join(random_dir, p))
    return [Ifg(join(random_dir, p)) for p in bases], random_dir


def same_exts_ifgs():
    """Return pair of Ifgs with same extents"""
    return [Ifg(join(PREP_TEST_TIF, f)) for f in ('0.tif', '1.tif')]


def test_extents_from_params():
    xf, yf = 1.0, 2.0
    xl, yl = 5.0, 7.0
    pars = {cf.IFG_XFIRST: xf, cf.IFG_XLAST: xl,
            cf.IFG_YFIRST: yf, cf.IFG_YLAST: yl}

    assert extents_from_params(pars) == CustomExts(xf, yf, xl, yl)


class PrepifgOutputTests(unittest.TestCase):
    """Tests aspects of the prepifg.py script, such as resampling."""

    def __init__(self, *args, **kwargs):
        super(PrepifgOutputTests, self).__init__(*args, **kwargs)

    @staticmethod
    def assert_geotransform_equal(files):
        """
        Asserts geotransforms for the given files are equivalent. Files can be paths
        to datasets, or GDAL dataset objects.
        """
        assert len(files) > 1, "Need more than 1 file to compare"
        if not all([hasattr(f, "GetGeoTransform") for f in files]):
            datasets = [gdal.Open(f) for f in files]
            assert all(datasets)
        else:
            datasets = files

        transforms = [ds.GetGeoTransform() for ds in datasets]
        head = transforms[0]
        for t in transforms[1:]:
            assert_array_almost_equal(t, head, decimal=6,
                                      err_msg="Extents do not match!")

    def setUp(self):
        self.xs = 0.000833333
        self.ys = -self.xs
        self.ifgs, self.random_dir = diff_exts_ifgs()
        self.ifg_paths = [i.data_path for i in self.ifgs]
        paths = ["geo_060619-061002_1rlks_1cr.tif",
                 "geo_060619-061002_1rlks_2cr.tif",
                 "geo_060619-061002_1rlks_3cr.tif",
                 "geo_060619-061002_4rlks_3cr.tif",
                 "geo_070326-070917_1rlks_1cr.tif",
                 "geo_070326-070917_1rlks_2cr.tif",
                 "geo_070326-070917_1rlks_3cr.tif",
                 "geo_070326-070917_4rlks_3cr.tif"]
        self.exp_files = [join(self.random_dir, p) for p in paths]

    @staticmethod
    def test_mlooked_paths():
        test_mlooked_path()

    @staticmethod
    def test_extents_from_params():
        test_extents_from_params()

    def tearDown(self):
        for exp_file in self.exp_files:
            if exists(exp_file):
                os.remove(exp_file)
        for ifg in self.ifgs:
            ifg.close()
        shutil.rmtree(self.random_dir)

    def _custom_ext_latlons(self):
        return [150.91 + (7 * self.xs),  # xfirst
                -34.17 + (16 * self.ys),  # yfirst
                150.91 + (27 * self.xs),  # 20 cells from xfirst
                -34.17 + (44 * self.ys)]  # 28 cells from yfirst

    def _custom_extents_tuple(self):
        return CustomExts(*self._custom_ext_latlons())

    def assert_projection_equal(self, files):
        """
        Asserts preojections for the given files are equivalent.
        Files can be paths to datasets, or GDAL dataset objects.
        """
        assert len(files) > 1, "Need more than 1 file to compare"
        if not all([hasattr(f, "GetGeoTransform") for f in files]):
            datasets = [gdal.Open(f) for f in files]
            assert all(datasets)
        else:
            datasets = files

        projections = [ds.GetProjection() for ds in datasets]
        head = projections[0]
        for t in projections[1:]:
            self.assertEqual(t, head)

    def test_multilooked_projection_same_as_geotiff(self):
        xlooks = ylooks = 1
        prepare_ifgs(self.ifg_paths, MAXIMUM_CROP, xlooks, ylooks)
        mlooked_paths = [mlooked_path(f, crop_out=MAXIMUM_CROP, looks=xlooks)
                         for f in self.ifg_paths]
        self.assert_projection_equal(self.ifg_paths + mlooked_paths)

    def test_default_max_extents(self):
        """Test ifgcropopt=2 crops datasets to max bounding box extents."""
        xlooks = ylooks = 1
        prepare_ifgs(self.ifg_paths, MAXIMUM_CROP, xlooks, ylooks)
        for f in [self.exp_files[1], self.exp_files[5]]:
            self.assertTrue(exists(f), msg="Output files not created")

        # output files should have same extents
        # NB: also verifies gdalwarp correctly copies geotransform across
        ifg = Ifg(self.exp_files[1])
        ifg.open()
        gt = ifg.dataset.GetGeoTransform()

        # copied from gdalinfo output
        exp_gt = (150.91, 0.000833333, 0, -34.17, 0, -0.000833333)
        for i, j in zip(gt, exp_gt):
            self.assertAlmostEqual(i, j)

        self.assert_geotransform_equal([self.exp_files[1], self.exp_files[5]])

        ifg.close()
        for i in self.ifgs:
            i.close()

    def test_min_extents(self):
        """Test ifgcropopt=1 crops datasets to min extents."""
        xlooks = ylooks = 1
        prepare_ifgs(self.ifg_paths, MINIMUM_CROP, xlooks, ylooks)
        ifg = Ifg(self.exp_files[0])
        ifg.open()

        # output files should have same extents
        # NB: also verifies gdalwarp correctly copies geotransform across
        # NB: expected data copied from gdalinfo output
        gt = ifg.dataset.GetGeoTransform()
        exp_gt = (150.911666666, 0.000833333, 0,
                  -34.172499999, 0, -0.000833333)
        for i, j in zip(gt, exp_gt):
            self.assertAlmostEqual(i, j)
        self.assert_geotransform_equal([self.exp_files[0], self.exp_files[4]])

        ifg.close()
        for i in self.ifgs:
            i.close()

    def test_custom_extents(self):
        xlooks = ylooks = 1
        cext = self._custom_extents_tuple()
        prepare_ifgs(self.ifg_paths, CUSTOM_CROP, xlooks, ylooks,
                     user_exts=cext)

        ifg = Ifg(self.exp_files[2])
        ifg.open()

        gt = ifg.dataset.GetGeoTransform()
        exp_gt = (cext.xfirst, self.xs, 0, cext.yfirst, 0, self.ys)

        for i, j in zip(gt, exp_gt):
            self.assertAlmostEqual(i, j)
        self.assert_geotransform_equal([self.exp_files[2], self.exp_files[6]])

        # close ifgs
        ifg.close()
        for i in self.ifgs:
            i.close()

    def test_exception_without_all_4_crop_parameters(self):
        """Test misaligned cropping extents raise errors."""
        xlooks = ylooks = 1
        # empty string and none raises exceptio
        for i in [None, '']:
            cext = (150.92, -34.18, 150.94, i)
            self.assertRaises(PreprocessError, prepare_ifgs, self.ifg_paths,
                              CUSTOM_CROP, xlooks, ylooks, user_exts=cext)
        # three parameters provided
        self.assertRaises(PreprocessError, prepare_ifgs, self.ifg_paths,
                          CUSTOM_CROP, xlooks, ylooks,
                          user_exts=(150.92, -34.18, 150.94))
        # close ifgs
        for i in self.ifgs:
            i.close()

    def test_custom_extents_misalignment(self):
        """Test misaligned cropping extents raise errors."""
        xlooks = ylooks = 1
        latlons = tuple(self._custom_ext_latlons())
        for i, _ in enumerate(['xfirst', 'yfirst', 'xlast', 'ylast']):
            # error = step / pi * [1000 100]
            for error in [0.265258, 0.026526]:
                tmp_latlon = list(latlons)
                tmp_latlon[i] += error
                cext = CustomExts(*tmp_latlon)

                self.assertRaises(PreprocessError, prepare_ifgs,
                                  self.ifg_paths, CUSTOM_CROP,
                                  xlooks, ylooks, user_exts=cext)
        # close ifgs
        for i in self.ifgs:
            i.close()

    def test_nodata(self):
        """Verify NODATA value copied correctly (amplitude band not copied)"""
        xlooks = ylooks = 1
        prepare_ifgs(self.ifg_paths, MINIMUM_CROP, xlooks, ylooks)

        for ex in [self.exp_files[0], self.exp_files[4]]:
            ifg = Ifg(ex)
            ifg.open()
            # NB: amplitude band doesn't have a NODATA value
            self.assertTrue(
                isnan(ifg.dataset.GetRasterBand(1).GetNoDataValue()))
            ifg.close()
        for i in self.ifgs:
            i.close()

    def test_nans(self):
        """Verify NaNs replace 0 in the multilooked phase band"""
        xlooks = ylooks = 1
        prepare_ifgs(self.ifg_paths, MINIMUM_CROP, xlooks, ylooks)

        for ex in [self.exp_files[0], self.exp_files[4]]:
            ifg = Ifg(ex)
            ifg.open()

            phase = ifg.phase_band.ReadAsArray()
            self.assertFalse((phase == 0).any())
            self.assertTrue((isnan(phase)).any())
            ifg.close()

        self.assertAlmostEqual(nanmax(phase), 4.247, 3)  # copied from gdalinfo
        self.assertAlmostEqual(nanmin(phase), 0.009, 3)  # copied from gdalinfo
        for i in self.ifgs:
            i.close()

    def test_multilook(self):
        """Test resampling method using a scaling factor of 4"""
        scale = 4  # assumes square cells
        self.ifgs.append(DEM(SYD_TEST_DEM_TIF))
        self.ifg_paths = [i.data_path for i in self.ifgs]
        cext = self._custom_extents_tuple()
        xlooks = ylooks = scale
        prepare_ifgs(self.ifg_paths, CUSTOM_CROP, xlooks, ylooks,
                     thresh=1.0, user_exts=cext)

        for n, ipath in enumerate([self.exp_files[3], self.exp_files[7]]):
            i = Ifg(ipath)
            i.open()
            self.assertEqual(i.dataset.RasterXSize, 20 / scale)
            self.assertEqual(i.dataset.RasterYSize, 28 / scale)

            # verify resampling
            path = join(PREP_TEST_TIF, "%s.tif" % n)
            ds = gdal.Open(path)
            src_data = ds.GetRasterBand(2).ReadAsArray()
            exp_resample = multilooking(src_data, scale, scale, thresh=0)
            self.assertEqual(exp_resample.shape, (7, 5))
            assert_array_almost_equal(exp_resample, i.phase_band.ReadAsArray())
            ds = None
            i.close()
            os.remove(ipath)

        # verify DEM has been correctly processed
        # ignore output values as resampling has already been tested for phase
        exp_dem_path = join(SYD_TEST_DEM_DIR, 'sydney_trimmed_4rlks_3cr.tif')
        self.assertTrue(exists(exp_dem_path))
        orignal_dem = DEM(SYD_TEST_DEM_TIF)
        orignal_dem.open()
        dem_dtype = orignal_dem.dataset.GetRasterBand(1).DataType
        orignal_dem.close()
        dem = DEM(exp_dem_path)
        dem.open()

        # test multilooked dem is of the same datatype as the original dem tif
        self.assertEqual(dem_dtype, dem.dataset.GetRasterBand(1).DataType)

        self.assertEqual(dem.dataset.RasterXSize, 20 / scale)
        self.assertEqual(dem.dataset.RasterYSize, 28 / scale)
        data = dem.height_band.ReadAsArray()
        self.assertTrue(data.ptp() != 0)
        # close ifgs
        dem.close()
        for i in self.ifgs:
            i.close()
        os.remove(exp_dem_path)

    def test_output_datatype(self):
        """Test resampling method using a scaling factor of 4"""
        scale = 4  # assumes square cells
        self.ifgs.append(DEM(SYD_TEST_DEM_TIF))
        self.ifg_paths = [i.data_path for i in self.ifgs] + [SYD_TEST_DEM_TIF]
        cext = self._custom_extents_tuple()
        xlooks = ylooks = scale
        prepare_ifgs(self.ifg_paths, CUSTOM_CROP, xlooks, ylooks,
                     thresh=1.0, user_exts=cext)

        for i in self.ifg_paths:
            mlooked_ifg = mlooked_path(i, xlooks, CUSTOM_CROP)
            ds1 = DEM(mlooked_ifg)
            ds1.open()
            ds2 = DEM(i)
            ds2.open()
            self.assertEqual(ds1.dataset.GetRasterBand(1).DataType,
                             ds2.dataset.GetRasterBand(1).DataType)
            ds1 = ds2 = None

    def test_invalid_looks(self):
        """Verify only numeric values can be given for multilooking"""
        values = [0, -1, -10, -100000.6, ""]
        for v in values:
            self.assertRaises(PreprocessError, prepare_ifgs, self.ifg_paths,
                              CUSTOM_CROP, xlooks=v, ylooks=1)

            self.assertRaises(PreprocessError, prepare_ifgs, self.ifg_paths,
                              CUSTOM_CROP, xlooks=1, ylooks=v)


class ThresholdTests(unittest.TestCase):
    """Tests for threshold of data -> NaN during resampling."""

    def test_nan_threshold_inputs(self):
        data = ones((1, 1))
        for thresh in [-10, -1, -0.5, 1.000001, 10]:
            self.assertRaises(ValueError, resample, data, 2, 2, thresh)

    @staticmethod
    def test_nan_threshold():
        # test threshold based on number of NaNs per averaging tile
        data = ones((2, 10))
        data[0, 3:] = nan
        data[1, 7:] = nan

        # key: NaN threshold as a % of pixels, expected result
        expected = [(0.0, [1, nan, nan, nan, nan]),
                    (0.25, [1, nan, nan, nan, nan]),
                    (0.5, [1, 1, nan, nan, nan]),
                    (0.75, [1, 1, 1, nan, nan]),
                    (1.0, [1, 1, 1, 1, nan])]

        for thresh, exp in expected:
            res = resample(data, xscale=2, yscale=2, thresh=thresh)
            assert_array_equal(res, reshape(exp, res.shape))

    @staticmethod
    def test_nan_threshold_alt():
        # test threshold on odd numbers
        data = ones((3, 6))
        data[0] = nan
        data[1, 2:5] = nan

        expected = [(0.4, [nan, nan]), (0.5, [1, nan]), (0.7, [1, 1])]
        for thresh, exp in expected:
            res = resample(data, xscale=3, yscale=3, thresh=thresh)
            assert_array_equal(res, reshape(exp, res.shape))


class SameSizeTests(unittest.TestCase):
    """Tests aspects of the prepifg.py script, such as resampling."""

    def __init__(self, *args, **kwargs):
        super(SameSizeTests, self).__init__(*args, **kwargs)
        self.xs = 0.000833333
        self.ys = -self.xs

    # TODO: check output files for same extents?
    # TODO: make prepifg dir readonly to test output to temp dir
    # TODO: move to class for testing same size option?
    def test_already_same_size(self):
        # should do nothing as layers are same size & no multilooking required
        ifgs = same_exts_ifgs()
        ifg_data_paths = [d.data_path for d in ifgs]
        res_tup = prepare_ifgs(ifg_data_paths, ALREADY_SAME_SIZE, 1, 1)
        res = [r[1] for r in res_tup]
        self.assertTrue(all(res))

    def test_already_same_size_mismatch(self):
        ifgs, random_dir = diff_exts_ifgs()
        ifg_data_paths = [d.data_path for d in ifgs]
        self.assertRaises(PreprocessError, prepare_ifgs,
                          ifg_data_paths, ALREADY_SAME_SIZE, 1, 1)
        for i in ifgs:
            i.close()
        shutil.rmtree(random_dir)

    # TODO: ensure multilooked files written to output dir
    def test_same_size_multilooking(self):
        ifgs = same_exts_ifgs()
        ifg_data_paths = [d.data_path for d in ifgs]
        xlooks = ylooks = 2
        prepare_ifgs(ifg_data_paths, ALREADY_SAME_SIZE, xlooks, ylooks)

        looks_paths = [mlooked_path(d, looks=xlooks, crop_out=ALREADY_SAME_SIZE)
                       for d in ifg_data_paths]
        mlooked = [Ifg(i) for i in looks_paths]
        for m in mlooked:
            m.open()
        self.assertEqual(len(mlooked), 2)

        for ifg in mlooked:
            self.assertAlmostEqual(ifg.x_step, xlooks * self.xs)
            self.assertAlmostEqual(ifg.x_step, ylooks * self.xs)
            # os.remove(ifg.data_path)


def test_mlooked_path():
    path = 'geo_060619-061002.tif'
    assert mlooked_path(path, looks=2, crop_out=4) == \
        'geo_060619-061002_2rlks_4cr.tif'

    path = 'some/dir/geo_060619-061002.tif'
    assert mlooked_path(path, looks=4, crop_out=2) == \
        'some/dir/geo_060619-061002_4rlks_2cr.tif'

    path = 'some/dir/geo_060619-061002_4rlks.tif'
    assert mlooked_path(path, looks=4, crop_out=8) == \
        'some/dir/geo_060619-061002_4rlks_4rlks_8cr.tif'


# class LineOfSightTests(unittest.TestCase):
# def test_los_conversion(self):
# TODO: needs LOS matrix
# TODO: this needs to work from config and incidence files on disk
# TODO: is convflag (see 'ifgconv' setting) used or just defaulted?
# TODO: los conversion has 4 options: 1: ignore, 2: vertical, 3: N/S, 4: E/W
# also have a 5th option of arbitrary azimuth angle (Pirate doesn't have this)
#    params = _default_extents_param()
#    params[IFG_CROP_OPT] = MINIMUM_CROP
#    params[PROJECTION_FLAG] = None
#    prepare_ifgs(params)


# def test_phase_conversion(self):
# TODO: check output data is converted to mm from radians (in prepifg??)
# raise NotImplementedError


class LocalMultilookTests(unittest.TestCase):
    """Tests for local testing functions"""

    @staticmethod
    def test_multilooking_thresh():
        data = ones((3, 6))
        data[0] = nan
        data[1, 2:5] = nan
        expected = [(6, [nan, nan]),
                    (5, [1, nan]),
                    (4, [1, 1])]
        scale = 3
        for thresh, exp in expected:
            res = multilooking(data, scale, scale, thresh)
            assert_array_equal(res, reshape(exp, res.shape))


def multilooking(src, xscale, yscale, thresh=0):
    """
    Port of looks.m from MATLAB Pirate.

    src: numpy array of phase data
    thresh: min number of non-NaNs required for a valid tile resampling
    """
    thresh = int(thresh)
    num_cells = xscale * yscale
    if thresh > num_cells or thresh < 0:
        msg = "Invalid threshold: %s (need 0 <= thr <= %s" % (thresh, num_cells)
        raise ValueError(msg)

    rows, cols = src.shape
    rows_lowres = int(floor(rows / yscale))
    cols_lowres = int(floor(cols / xscale))
    dest = ones((rows_lowres, cols_lowres)) * nan

    size = xscale * yscale
    for row in range(rows_lowres):
        for col in range(cols_lowres):
            ys = row * yscale
            ye = ys + yscale
            xs = col * xscale
            xe = xs + xscale

            patch = src[ys:ye, xs:xe]
            num_values = num_cells - npsum(isnan(patch))

            if num_values >= thresh and num_values > 0:
                # nanmean() only works on 1g axis
                reshaped = patch.reshape(size)
                dest[row, col] = nanmean(reshaped)

    return dest


class MatlabEqualityTestRoipacSydneyTestData(unittest.TestCase):
    """
    Matlab to python roipac prepifg equality test for sydney test data
    """

    def setUp(self):
        from tests.common import sydney_data_setup
        self.ifgs = sydney_data_setup()
        self.ifg_paths = [i.data_path for i in self.ifgs]
        prepare_ifgs(self.ifg_paths, crop_opt=1, xlooks=1, ylooks=1)
        looks_paths = [mlooked_path(d, looks=1, crop_out=1)
                       for d in self.ifg_paths]
        self.ifgs_with_nan = [Ifg(i) for i in looks_paths]
        for ifg in self.ifgs_with_nan:
            ifg.open()

    def tearDown(self):
        for i in self.ifgs_with_nan:
            if os.path.exists(i.data_path):
                i.close()
                os.remove(i.data_path)

    def test_matlab_prepifg_equality_array(self):
        """
        Matlab to python prepifg equality test
        """
        # path to csv folders from matlab output
        onlyfiles = [
            fln for fln in os.listdir(SYD_TEST_MATLAB_PREPIFG_DIR)
            if os.path.isfile(os.path.join(SYD_TEST_MATLAB_PREPIFG_DIR, fln))
            and fln.endswith('.csv') and fln.__contains__('_rad_')
            ]

        for fln in onlyfiles:
            ifg_data = np.genfromtxt(os.path.join(
                SYD_TEST_MATLAB_PREPIFG_DIR, fln), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if fln.split('_rad_')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('.')[0]:
                    np.testing.assert_array_almost_equal(ifg_data,
                                                         self.ifgs_with_nan[
                                                             k].phase_data,
                                                         decimal=2)

    def test_matlab_prepifg_and_convert_phase(self):
        """
        Matlab to python prepifg equality test
        """
        # path to csv folders from matlab output
        for i in self.ifgs_with_nan:
            if not i.mm_converted:
                i.convert_to_mm()
        onlyfiles = [
            f for f in os.listdir(SYD_TEST_MATLAB_PREPIFG_DIR)
            if os.path.isfile(os.path.join(SYD_TEST_MATLAB_PREPIFG_DIR, f))
            and f.endswith('.csv') and f.__contains__('_mm_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                SYD_TEST_MATLAB_PREPIFG_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_mm_')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_unw.')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(
                        ifg_data, self.ifgs_with_nan[k].phase_data, decimal=2)

                    # means must also be equal
                    self.assertAlmostEqual(
                        np.nanmean(ifg_data),
                        np.nanmean(self.ifgs_with_nan[k].phase_data),
                        places=4)

                    # number of nans must equal
                    self.assertEqual(
                        np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs_with_nan[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))


class TestOneIncidenceOrElevationMap(unittest.TestCase):

    def setUp(self):
        self.base_dir = tempfile.mkdtemp()
        self.conf_file = tempfile.mktemp(suffix='.conf', dir=self.base_dir)
        self.ifgListFile = os.path.join(common.SYD_TEST_GAMMA, 'ifms_17')

    def tearDown(self):
        shutil.rmtree(self.base_dir)

    def make_input_files(self, inc='', ele=''):
        with open(self.conf_file, 'w') as conf:
            conf.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, common.SYD_TEST_GAMMA))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '1'))
            conf.write('{}: {}\n'.format(LUIGI, '0'))
            conf.write('{}: {}\n'.format(
                DEM_HEADER_FILE, os.path.join(
                    common.SYD_TEST_GAMMA, '20060619_utm_dem.par')))
            conf.write('{}: {}\n'.format(IFG_LKSX, '1'))
            conf.write('{}: {}\n'.format(IFG_LKSY, '1'))
            conf.write('{}: {}\n'.format(IFG_CROP_OPT, '1'))
            conf.write('{}: {}\n'.format(NO_DATA_AVERAGING_THRESHOLD, '0.5'))
            conf.write('{}: {}\n'.format(SLC_DIR, ''))
            conf.write('{}: {}\n'.format(DEM_FILE, common.SYD_TEST_DEM_GAMMA))
            conf.write('{}: {}\n'.format(APS_INCIDENCE_MAP, inc))
            conf.write('{}: {}\n'.format(APS_ELEVATION_MAP, ele))
            conf.write('{}: {}\n'.format(APS_CORRECTION, '1'))
            conf.write('{}: {}\n'.format(APS_METHOD, '2'))

    def test_only_inc_file_created(self):
        inc_ext = 'inc'
        ele_ext = 'lv_theta'
        self.make_input_files(inc=common.SYD_TEST_INCIDENCE)
        self.common_check(inc_ext, ele_ext)

    def test_only_ele_file_created(self):
        inc_ext = 'inc'
        ele_ext = 'lv_theta'
        self.make_input_files(ele=common.SYD_TEST_ELEVATION)
        self.common_check(ele_ext, inc_ext)

    def common_check(self, ele, inc):
        os.path.exists(self.conf_file)
        params = cf.get_config_params(self.conf_file)
        sys.argv = ['dummy', self.conf_file]
        run_prepifg.main(params)
        # test 17 geotiffs created
        geotifs = glob.glob(os.path.join(self.base_dir, '*_unw.tif'))
        self.assertEqual(17, len(geotifs))
        # test dem geotiff created
        demtif = glob.glob(os.path.join(self.base_dir, '*_dem.tif'))
        self.assertEqual(1, len(demtif))
        # elevation/incidence file
        ele = glob.glob(os.path.join(self.base_dir,
                                     '*utm_{ele}.tif'.format(ele=ele)))[0]
        self.assertTrue(os.path.exists(ele))
        # mlooked tifs
        mlooked_tifs = [f for f in
                        glob.glob(os.path.join(self.base_dir, '*.tif'))
                        if "cr" in f and "rlks" in f]
        # 19 including 17 ifgs, 1 dem and one incidence
        self.assertEqual(19, len(mlooked_tifs))
        inc = glob.glob(os.path.join(self.base_dir,
                                     '*utm_{inc}.tif'.format(inc=inc)))
        self.assertEqual(0, len(inc))


if __name__ == "__main__":
    unittest.main()
