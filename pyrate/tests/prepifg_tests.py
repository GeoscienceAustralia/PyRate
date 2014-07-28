'''
Tests for prepifg.py: resampling, subsetting etc

Ben Davies, NCI
'''

import os, sys
import unittest
from os.path import exists, join
from collections import namedtuple

from math import floor
from scipy.stats.stats import nanmean
from numpy import isnan, nanmax, nanmin
from numpy import ones, nan, reshape, sum as npsum
from numpy.testing import assert_array_almost_equal, assert_array_equal

from pyrate.prepifg import CUSTOM_CROP, MAXIMUM_CROP, MINIMUM_CROP, ALREADY_SAME_SIZE
from pyrate.prepifg import prepare_ifgs, resample, PreprocessError

from pyrate.shared import Ifg, DEM
from pyrate.config import OBS_DIR
from pyrate.config import DEM_FILE, IFG_FILE_LIST

from common import PREP_TEST_TIF, SYD_TEST_DEM_DIR, SYD_TEST_DEM_TIF

import gdal
gdal.UseExceptions()

if not exists(PREP_TEST_TIF):
	sys.exit("ERROR: Missing 'prepifg' dir for unittests\n")


CustomExts = namedtuple('CustExtents', ['xfirst', 'yfirst', 'xlast', 'ylast'])



def _default_extents_param():
	# create dummy params file (relative paths to prevent chdir calls)
	return { IFG_FILE_LIST: join(PREP_TEST_TIF, 'ifms'),
			OBS_DIR: PREP_TEST_TIF }


class PrepifgOutputTests(unittest.TestCase):
	"""Tests aspects of the prepifg.py script, such as resampling."""

	def __init__(self, *args, **kwargs):
		super(PrepifgOutputTests, self).__init__(*args, **kwargs)

	def setUp(self):
		self.xs = 0.000833333
		self.ys = -self.xs
		paths = ["geo_060619-061002_1rlks.tif", "geo_070326-070917_1rlks.tif"]
		self.exp_files = [join(PREP_TEST_TIF, p) for p in paths]

	def tearDown(self):
		# clear temp output files after each run
		for f in self.exp_files:
			if exists(f):
				os.remove(f)


	def _custom_extents_param(self):
		'Convenience function to create custom cropping extents params'
		params = {}
		params[IFG_FILE_LIST] = join(PREP_TEST_TIF, 'ifms')
		params[OBS_DIR] = PREP_TEST_TIF
		return params

	def _cust_ext_latlons(self):
		return [150.91 + (7 * self.xs), # xfirst
				-34.17 + (16 * self.ys), # yfirst
				150.91 + (27 * self.xs), # 20 cells from xfirst
				-34.17 + (44 * self.ys)] # 28 cells from yfirst

	def _custom_extents_tuple(self):
		return CustomExts(*self._cust_ext_latlons())

	def test_default_max_extents(self):
		'Test ifgcropopt=2 gives datasets cropped to max extents bounding box.'
		params = _default_extents_param()
		xlooks = ylooks = 1
		prepare_ifgs(MAXIMUM_CROP, xlooks, ylooks, params)
		for f in self.exp_files:
			self.assertTrue(exists(f), msg="Output files not created")

		# output files should have same extents
		# NB: also verifies gdalwarop correctly copies geotransform across
		ifg = Ifg(self.exp_files[0])
		ifg.open()
		gt = ifg.dataset.GetGeoTransform()

		# copied from gdalinfo output
		exp_gt = (150.91, 0.000833333, 0, -34.17, 0, -0.000833333)
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)


	def test_min_extents(self):
		"""Test ifgcropopt=1 crops datasets to min extents."""
		params = _default_extents_param()
		xlooks = ylooks = 1
		prepare_ifgs(MINIMUM_CROP, xlooks, ylooks, params)
		ifg = Ifg(self.exp_files[0])
		ifg.open()

		# output files should have same extents
		# NB: also verifies gdalwarp correctly copies geotransform across
		# NB: expected data copied from gdalinfo output
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (150.911666666, 0.000833333, 0, -34.172499999, 0, -0.000833333)
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)


	def test_custom_extents(self):
		params = self._custom_extents_param()
		cext = self._custom_extents_tuple()
		xlooks = ylooks = 1

		prepare_ifgs(CUSTOM_CROP, xlooks, ylooks, params, custom_crop=cext)

		ifg = Ifg(self.exp_files[0])
		ifg.open()

		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (cext.xfirst, self.xs, 0, cext.yfirst, 0, self.ys)

		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)

		assert_geotransform_equal(self.exp_files)


	def test_custom_extents_misalignment(self):
		"""Test misaligned cropping extents raise errors."""
		params = self._custom_extents_param()
		xlooks = ylooks = 1

		latlons = tuple(self._cust_ext_latlons())
		for i, _ in enumerate(['xfirst', 'yfirst', 'xlast', 'ylast']):
			for error in [0.1, 0.001, 0.0001, 0.00001, 0.000001]:
				tmp_latlon = list(latlons)
				tmp_latlon[i] += error
				cext = CustomExts(*tmp_latlon)

				self.assertRaises(PreprocessError, prepare_ifgs, CUSTOM_CROP,
								xlooks, ylooks, params, custom_crop=cext)

	def test_nodata(self):
		'Verify NODATA value is copied correctly (amplitude band not copied)'
		params = _default_extents_param()
		xlooks = ylooks = 1
		prepare_ifgs(MINIMUM_CROP, xlooks, ylooks, params)

		for ex in self.exp_files:
			ifg = Ifg(ex)
			ifg.open()
			# NB: amplitude band doesn't have a NODATA value
			self.assertTrue(isnan(ifg.dataset.GetRasterBand(1).GetNoDataValue()))

	def test_nans(self):
		"""Verify that NaNs replace 0 in the multilooked phase band"""
		params = _default_extents_param()
		xlooks = ylooks = 1
		prepare_ifgs(MINIMUM_CROP, xlooks, ylooks, params)

		for ex in self.exp_files:
			ifg = Ifg(ex)
			ifg.open()

			phase = ifg.phase_band.ReadAsArray()
			self.assertFalse((phase == 0).any() )
			self.assertTrue((isnan(phase)).any() )

		self.assertAlmostEqual(nanmax(phase), 4.247, 3) # copied from gdalinfo
		self.assertAlmostEqual(nanmin(phase), 0.009, 3) # copied from gdalinfo


	def test_multilook(self):
		"""Test resampling method using a scaling factor of 4"""
		scale = 4 # assumes square cells
		params = self._custom_extents_param()
		cext = self._custom_extents_tuple()

		xlooks = scale
		ylooks = scale
		params[DEM_FILE] = SYD_TEST_DEM_TIF
		prepare_ifgs(CUSTOM_CROP, xlooks, ylooks, params,
					thresh=1.0, custom_crop=cext) # if all nans, ignore cell

		# check file names have been updated
		for f in self.exp_files:
			self.assertFalse(exists(f))

		exp_files = [s.replace('_1r', '_%sr' % scale) for s in self.exp_files]

		for n, ipath in enumerate(exp_files):
			i = Ifg(ipath)
			i.open()
			self.assertEqual(i.dataset.RasterXSize, 20 / scale)
			self.assertEqual(i.dataset.RasterYSize, 28 / scale)

			# verify resampling
			path = join(PREP_TEST_TIF, "%s.tif" % n)
			ds = gdal.Open(path)
			src_data = ds.GetRasterBand(2).ReadAsArray()
			exp_resample = multilooking(src_data, scale, scale, thresh=0)
			self.assertEqual(exp_resample.shape, (7,5))
			assert_array_almost_equal(exp_resample, i.phase_band.ReadAsArray())
			os.remove(ipath)

		# verify DEM has been correctly processed
		# ignore output values as resampling has already been tested for phase
		exp_dem_path = join(SYD_TEST_DEM_DIR, 'sydney_trimmed_4rlks.dem')
		self.assertTrue(exists(exp_dem_path))

		dem = DEM(exp_dem_path)
		dem.open()
		self.assertEqual(dem.dataset.RasterXSize, 20 / scale)
		self.assertEqual(dem.dataset.RasterYSize, 28 / scale)
		data = dem.height_band.ReadAsArray()
		self.assertTrue(data.ptp() != 0)
		os.remove(dem.data_path)


	def test_invalid_looks(self):
		"""Verify only numeric values can be given for multilooking"""
		params = self._custom_extents_param()
		values = [0, -1, -10, -100000.6, ""]
		for v in values:
			self.assertRaises(PreprocessError, prepare_ifgs, CUSTOM_CROP,
								xlooks=v, ylooks=1, params=params)

			self.assertRaises(PreprocessError, prepare_ifgs, CUSTOM_CROP,
								xlooks=1, ylooks=v, params=params)


	def test_nan_threshold_inputs(self):
		data = ones((1,1))
		for thresh in [-10, -1, -0.5, 1.000001, 10]:
			self.assertRaises(ValueError, resample, data, 2, 2, thresh)


	def test_nan_threshold(self):
		# test threshold based on number of NaNs per averaging tile
		data = ones((2,10))
		data[0,3:] = nan
		data[1,7:] = nan

		# key: NaN threshold as a % of pixels, expected result
		expected = [ (0.0, [1, nan, nan, nan, nan]),
					(0.25, [1, nan, nan, nan, nan]),
					(0.5, [1, 1, nan, nan, nan]),
					(0.75, [1, 1, 1, nan, nan]),
					(1.0, [1, 1, 1, 1, nan])  ]

		for thresh, exp in expected:
			res = resample(data, xscale=2, yscale=2, thresh=thresh)
			assert_array_equal(res, reshape(exp, res.shape))


	def test_nan_threshold_alt(self):
		# test threshold on odd numbers
		data = ones((3,6))
		data[0] = nan
		data[1,2:5] = nan

		expected = [ (0.4, [nan,nan]), (0.5, [1, nan]), (0.7, [1, 1]) ]
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
		params = _default_extents_param()
		params[IFG_FILE_LIST] = join(PREP_TEST_TIF, 'ifms2')
		self.assertEqual(prepare_ifgs(ALREADY_SAME_SIZE, 1, 1, params), None)

	def test_already_same_size_mismatch(self):
		params = _default_extents_param()
		self.assertRaises(PreprocessError, prepare_ifgs, ALREADY_SAME_SIZE, 1, 1, params)

	# TODO: ensure multilooked files written to output dir
	def test_same_size_multilooking(self):
		params = _default_extents_param()
		params[IFG_FILE_LIST] = join(PREP_TEST_TIF, 'ifms2')
		xlooks = 2
		ylooks = 2

		mlooked = prepare_ifgs(ALREADY_SAME_SIZE, xlooks, ylooks, params)
		self.assertEqual(len(mlooked), 2)

		for ifg in mlooked:
			self.assertEqual(ifg.x_step, xlooks * self.xs)
			self.assertEqual(ifg.x_step, ylooks * self.xs)
			os.remove(ifg.data_path)


#class LineOfSightTests(unittest.TestCase):
	#def test_los_conversion(self):
		# TODO: needs LOS matrix
		# TODO: this needs to work from config and incidence files on disk
		# TODO: is convflag (see 'ifgconv' setting) used or just defaulted?
		# TODO: los conversion has 4 options: 1: ignore, 2: vertical, 3: N/S, 4: E/W
		# also have a 5th option of arbitrary azimuth angle (Pirate doesn't have this)
	#	params = _default_extents_param()
	#	params[IFG_CROP_OPT] = MINIMUM_CROP
	#	params[PROJECTION_FLAG] = None
	#	prepare_ifgs(params)


	#def test_phase_conversion(self):
		# TODO: check output data is converted to mm from radians (in prepifg??)
		#raise NotImplementedError



class LocalMultilookTests(unittest.TestCase):
	'''Tests for local testing functions'''

	def test_multilooking_thresh(self):
		data = ones((3,6))
		data[0] = nan
		data[1,2:5] = nan
		expected = [ (6, [nan,nan]),
					(5, [1, nan]),
					(4, [1, 1]) ]
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
	for r in range(rows_lowres):
		for c in range(cols_lowres):
			ys = r * yscale
			ye = ys + yscale
			xs = c * xscale
			xe = xs + xscale

			patch = src[ys:ye, xs:xe]
			num_values = num_cells - npsum(isnan(patch))

			if num_values >= thresh and num_values > 0:
				reshaped = patch.reshape(size) # NB: nanmean only works on 1 axis
				dest[r,c] = nanmean(reshaped)

	return dest


def assert_geotransform_equal(files):
	"""
	Asserts geotransforms for the given files are equivalent. Files can be paths
	to datasets, or GDAL dataset objects.
	"""
	assert len(files) > 1, "Need more than 1 file to compare"
	if not all( [hasattr(f, "GetGeoTransform") for f in files] ):
		datasets = [gdal.Open(f) for f in files]
		assert all(datasets)
	else:
		datasets = files

	transforms = [ds.GetGeoTransform() for ds in datasets]
	head = transforms[0]
	for t in transforms[1:]:
		assert t == head, "Extents do not match!"
