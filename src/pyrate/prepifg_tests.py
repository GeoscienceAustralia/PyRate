# Tests for prepifg.py
#


import os
import unittest
from os.path import exists, join

try:
	from osgeo import gdal
except:
	import gdal

import prepifg
from shared import Ifg
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST
from config import IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST


class OutputTests(unittest.TestCase):
	"""TODO"""
	
	def setUp(self):
		self.testdir = "../../tests/prepifg"
		tmp = ["obs/geo_060619-061002.unw.rsc", "obs/geo_070326-070917.unw.rsc"]
		self.hdr_files = [join(self.testdir, t) for t in tmp]
		
		tmp = ["obs/geo_060619-061002_TODO.tif", "obs/geo_070326-070917_TODO.tif"]	
		self.exp_files = [join(self.testdir, t) for t in tmp]
		
		# clear tmp test files before each run
		for f in self.exp_files:
			if exists(f):
				os.remove(f)


	def test_default_max_extents(self):
		"""Test ifgcropopt=2 gives datasets cropped to max extents bounding box."""
		
		# create dummy params file (relative paths to prevent chdir calls)
		params = {IFG_CROP_OPT: prepifg.MAXIMUM_CROP, IFG_LKSX: 1, IFG_LKSY: 1}
		params[IFG_FILE_LIST] = join(self.testdir, 'obs/ifms')
		params[OBS_DIR] = join(self.testdir,"obs/")
		
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
		
		# create dummy params file (relative paths to prevent chdir calls)
		params = {IFG_CROP_OPT: prepifg.MINIMUM_CROP, IFG_LKSX: 1, IFG_LKSY: 1}
		params[IFG_FILE_LIST] = join(self.testdir, 'obs/ifms')
		params[OBS_DIR] = join(self.testdir,"obs/")
		
		# output files should have same extents
		prepifg.prepare_ifgs(params)
		ifg = Ifg(self.exp_files[0], self.hdr_files[0])
		ifg.open()
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (150.911666666, 0.000833333, 0, -34.174166665, 0, -0.000833333) # copied from gdalinfo output
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)
	
	
	def test_custom_extents(self):
		# create dummy params file (relative paths to prevent chdir calls)
		params = {IFG_CROP_OPT: prepifg.CUSTOM_CROP, IFG_LKSX: 1, IFG_LKSY: 1}
		xs = 0.000833333
		ys = -xs
		params[IFG_XFIRST] = 150.91 + (7 * xs)
		params[IFG_YFIRST] = -34.17 + (16 * ys)
		params[IFG_XLAST] = 150.9491667 + xs
		params[IFG_YLAST] = -34.23 - (4 * xs)
		params[IFG_FILE_LIST] = join(self.testdir, 'obs/ifms')
		params[OBS_DIR] = join(self.testdir,"obs/")
		
		# output files should have same extents
		prepifg.prepare_ifgs(params)
		ifg = Ifg(self.exp_files[0], self.hdr_files[0])
		ifg.open()
		gt = ifg.dataset.GetGeoTransform()
		exp_gt = (params[IFG_XFIRST], xs, 0, params[IFG_YFIRST], 0, ys)
		for i,j in zip(gt, exp_gt):
			self.assertAlmostEqual(i, j)
		assert_geotransform_equal(self.exp_files)		


	def test_nodata(self):
		# TODO: ensure the nodata values are copied correctly for each band
		pass



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
