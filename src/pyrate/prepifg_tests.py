# Tests for prepifg.py
#


import os
import unittest
from os.path import exists, join

from numpy import  isnan, nanmax, nanmin

try:
	from osgeo import gdal
except:
	import gdal

gdal.UseExceptions()

import prepifg
from shared import Ifg
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST
from config import IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST


class OutputTests(unittest.TestCase):
	"""TODO"""
	
	def setUp(self):
		self.xs = 0.000833333
		self.ys = -self.xs
		
		self.testdir = "../../tests/prepifg"
		tmp = ["obs/geo_060619-061002.unw.rsc", "obs/geo_070326-070917.unw.rsc"]
		self.hdr_files = [join(self.testdir, t) for t in tmp]
		
		tmp = ["obs/geo_060619-061002_TODO.tif", "obs/geo_070326-070917_TODO.tif"]	
		self.exp_files = [join(self.testdir, t) for t in tmp]


	def tearDown(self):
		# clear tmp output files after each run
		for f in self.exp_files:
			if exists(f):
				os.remove(f)


	def _custom_extents_param(self):
		"""Convenience function to create custom cropping extents params"""
		
		params = {IFG_CROP_OPT: prepifg.CUSTOM_CROP, IFG_LKSX: 1, IFG_LKSY: 1}		
		params[IFG_XFIRST] = 150.91 + (7 * self.xs)
		params[IFG_YFIRST] = -34.17 + (16 * self.ys)
		params[IFG_XLAST] = 150.91 + (27 * self.xs) # 20 cells across from X_FIRST
		params[IFG_YLAST] = -34.17 + (44 * self.ys) # 28 cells across from Y_FIRST
		params[IFG_FILE_LIST] = join(self.testdir, 'obs/ifms')
		params[OBS_DIR] = join(self.testdir,"obs/")
		return params


	def _default_extents_param(self):
		# create dummy params file (relative paths to prevent chdir calls)
		return {IFG_LKSX: 1, IFG_LKSY: 1, IFG_FILE_LIST: join(self.testdir, 'obs/ifms'),
						OBS_DIR: join(self.testdir,"obs/") }	


	def test_default_max_extents(self):
		"""Test ifgcropopt=2 gives datasets cropped to max extents bounding box."""

		params = self._default_extents_param()
		params[IFG_CROP_OPT] = prepifg.MAXIMUM_CROP 
		prepifg.prepare_ifgs(params)
		for f in self.exp_files:
			self.assertTrue(exists(f), msg="Output files not created")

		# output files should have same extents
		ifg = Ifg(self.exp_files[0], self.hdr_files[0])
		ifg.open()
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (150.91, 0.000833333, 0, -34.17, 0, -0.000833333) # copied from gdalinfo output
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)


	def test_min_extents(self):
		"""Test ifgcropopt=1 crops datasets to min extents."""
		
		params = self._default_extents_param()		
		params[IFG_CROP_OPT] = prepifg.MINIMUM_CROP
		prepifg.prepare_ifgs(params)
		ifg = Ifg(self.exp_files[0], self.hdr_files[0])
		ifg.open()
		
		# output files should have same extents
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (150.911666666, 0.000833333, 0, -34.172499999, 0, -0.000833333) # copied from gdalinfo output
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)
	
	
	def test_custom_extents(self):
		params = self._custom_extents_param()
		prepifg.prepare_ifgs(params)
		ifg = Ifg(self.exp_files[0], self.hdr_files[0])
		ifg.open()
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (params[IFG_XFIRST], self.xs, 0, params[IFG_YFIRST], 0, self.ys)
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)
		
		# TODO: ensure extents are as specified?
		
	
	def test_custom_extents_misalignment(self):
		"""Test misaligned cropping extents raise errors."""
		
		for key in [IFG_XFIRST, IFG_YFIRST, IFG_XLAST, IFG_YLAST]:
			params = self._custom_extents_param() # reset params to prevent old params causing failure
			backup = params[key]
			
			# try different errors for each var
			for error in [0.1, 0.001, 0.0001, 0.00001, 0.000001]:
				params[key] = backup + error
				self.assertRaises(prepifg.PreprocessingException, prepifg.prepare_ifgs, params, use_exceptions=True)


	def test_nodata(self):
		"""Verify NODATA values are copied correctly for both bands"""
		params = self._default_extents_param()		
		params[IFG_CROP_OPT] = prepifg.MINIMUM_CROP
		prepifg.prepare_ifgs(params)
		
		for ex, hdr in zip(self.exp_files, self.hdr_files):
			ifg = Ifg(ex, hdr)
			ifg.open()	
			# NB: amplitude band doesn't have a NODATA value
			self.assertTrue(isnan(ifg.dataset.GetRasterBand(2).GetNoDataValue()))


	def test_nans(self):
		"""Verify that NaNs replace 0 in the multilooked phase band"""
		params = self._default_extents_param()		
		params[IFG_CROP_OPT] = prepifg.MINIMUM_CROP
		prepifg.prepare_ifgs(params)
		
		for ex, hdr in zip(self.exp_files, self.hdr_files):
			ifg = Ifg(ex, hdr)
			ifg.open()
			
			phase = ifg.phase_band.ReadAsArray()
			self.assertFalse((phase == 0).any() )
			self.assertTrue((isnan(phase)).any() )
		
		self.assertAlmostEqual(nanmax(phase), 4.247, 3) # copied from gdalinfo
		self.assertAlmostEqual(nanmin(phase), 0.009, 3) # copied from gdalinfo



def assert_geotransform_equal(files):
	"""Asserts geotransforms for the given files are equivalent. Files can be paths
	to datasets, or GDAL dataset objects."""
	
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
