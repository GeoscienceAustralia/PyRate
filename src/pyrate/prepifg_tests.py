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
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST



class OutputTests(unittest.TestCase):
	"""TODO"""
	
	def setUp(self):
		self.testdir = "../../tests/prepifg"
		tmp = ["obs/geo_060619-061002_TODO.tif", "obs/geo_070326-070917_TODO.tif"]	
		self.exp_files = [join(self.testdir, t) for t in tmp]
		
		# clear tmp test files before each run
		for f in self.exp_files:
			if exists(f):
				os.remove(f)
	

	def test_default_max_extents(self):
		"""Test ifgcropopt=2 gives datasets cropped to max extents bounding box."""
		
		# create dummy params file (large paths to prevent chdir calls)
		params = {IFG_CROP_OPT: 2, IFG_LKSX: 1, IFG_LKSY: 1}
		params[IFG_FILE_LIST] = join(self.testdir, 'obs/ifms')
		params[OBS_DIR] = join(self.testdir,"obs/")
		
		prepifg.prepare_ifgs(params)
		
		for f in self.exp_files:
			self.assertTrue(exists(f), msg="Output files not created")

		# should have 2 output files of the same bounds
		assert_geotransform_equal(self.exp_files)


	def test_nodata(self):
		# TODO: ensure the nodata values are copied correctly for each band
		pass



def assert_geotransform_equal(files):
	"""Asserts geotransforms for the given files are equivalent. Files can be paths
	to datasets, or GDAL dataset objects."""

	if not all( [hasattr(f, "GetGeoTransform") for f in files] ):
		datasets = [gdal.Open(f) for f in files]
		assert all(datasets)
	else:
		datasets = files

	transforms = [ds.GetGeoTransform() for ds in datasets]
	head = transforms[0]
	for t in transforms[1:]:
		assert t == head, "Extents do not match!" 
