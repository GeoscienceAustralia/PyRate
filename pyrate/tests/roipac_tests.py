'''
Tests for ROIPAC header translation module.

Created on 12/09/2012
@author: Ben Davies, NCI
         ben.davies@anu.edu.au
'''

import os
import sys
import unittest
from datetime import date
from os.path import exists, join

from pyrate import roipac
import pyrate.ifgconstants as ifc
from common import HEADERS_TEST_DIR, PREP_TEST_OBS, PREP_TEST_TIF
from common import SYD_TEST_DEM_UNW, SYD_TEST_DEM_HDR, SYD_TEST_DEM_DIR, SYD_TEST_OBS

from numpy.testing import assert_array_almost_equal

import gdal
from pyrate.roipac import RoipacException
gdal.UseExceptions()


#  sanity checking
if not exists(HEADERS_TEST_DIR):
	sys.exit("ERROR: Missing the 'headers' data for unittests\n")

# constants
SHORT_HEADER_PATH = join(SYD_TEST_OBS, 'geo_060619-061002.unw.rsc')
FULL_HEADER_PATH  = join(HEADERS_TEST_DIR, "geo_060619-060828.unw.rsc")



class RoipacToGeoTiffTests(unittest.TestCase):
	'Tests conversion of GAMMA rasters to custom PyRate GeoTIFF'
	
	@classmethod
	def setUpClass(cls):
		hdr_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw.rsc')
		cls.HDRS = roipac.parse_header(hdr_path)

	def test_to_geotiff_dem(self):
		hdr = roipac.parse_header(SYD_TEST_DEM_HDR)
		dest = "/tmp/tmp_roipac_dem.tif"

		roipac.to_geotiff(hdr, SYD_TEST_DEM_UNW, dest, nodata=0)
		exp_path = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.tif')
		exp_ds = gdal.Open(exp_path)
		ds = gdal.Open(dest)

		# compare data and geographic headers
		assert_array_almost_equal(exp_ds.ReadAsArray(), ds.ReadAsArray())
		self.compare_rasters(ds, exp_ds)
		self.assertIsNotNone(ds.GetMetadata())
		os.remove(dest)


	def test_to_geotiff_ifg(self):
		# tricker: needs ifg header, and DEM one for extents 
		hdrs = self.HDRS.copy()
		hdrs[ifc.PYRATE_DATUM] = 'WGS84'

		dest = '/tmp/tmp_roipac_ifg.tif'
		data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
		roipac.to_geotiff(hdrs, data_path, dest, nodata=0)		
				
		ds = gdal.Open(dest)
		band = ds.GetRasterBand(1)
		self.assertEqual(ds.RasterCount, 1)
		
		# TODO: result is only a single band dataset, make exp tif 1 band?
		exp_path = join(PREP_TEST_TIF, 'geo_060619-061002.tif')
		exp_ds = gdal.Open(exp_path)
		exp_band = exp_ds.GetRasterBand(1)
		
		# compare data and geographic headers
		assert_array_almost_equal(exp_band.ReadAsArray(), band.ReadAsArray())
		self.compare_rasters(ds, exp_ds)

		md = ds.GetMetadata()
		date1 = date(2006, 6, 19)
		date2 = date(2006, 10, 2)
		diff = (date2 - date1).days
		self.assertTrue(md[ifc.PYRATE_DATE] == str(date1))
		self.assertTrue(md[ifc.PYRATE_DATE2] == str(date2))
		self.assertTrue(md[ifc.PYRATE_TIME_SPAN] == str(diff / 365.25))
		
		wavelen = float(md[ifc.PYRATE_WAVELENGTH_METRES])
		self.assertAlmostEqual(wavelen, 0.0562356424)
		
	def test_to_geotiff_wrong_input_data(self):
		# ensure failure if TIF/other file used instead of binary UNW data
		dest = '/tmp/tmp_roipac_ifg.tif'
		data_path = join(PREP_TEST_TIF, 'geo_060619-061002.tif')
		self.assertRaises(roipac.RoipacException, roipac.to_geotiff,
							self.HDRS, data_path, dest, nodata=0)

	def test_bad_projection(self):
		hdrs = self.HDRS.copy()
		hdrs[ifc.PYRATE_DATUM] = 'bad datum string'
		dest = '/tmp/tmp_roipac_ifg2.tif'
		data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
		self.assertRaises(RoipacException, roipac.to_geotiff, hdrs, data_path, dest, 0)

	def test_mismatching_cell_resolution(self):
		hdrs = self.HDRS.copy()
		hdrs[ifc.PYRATE_X_STEP] = 0.1 # fake a mismatch
		hdrs[ifc.PYRATE_DATUM] = 'WGS84'
		data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
		dest = '/tmp/fake'
		
		self.assertRaises(RoipacException, roipac.to_geotiff, hdrs, data_path, dest, 0)

	def compare_rasters(self, ds, exp_ds):
		band = ds.GetRasterBand(1)
		exp_band = exp_ds.GetRasterBand(1)		

		nodata = band.GetNoDataValue()
		self.assertFalse(nodata is None)
		self.assertEqual(exp_band.GetNoDataValue(), nodata)
		
		pj = ds.GetProjection()		
		self.assertTrue('WGS 84' in pj)
		self.assertEqual(exp_ds.GetProjection(), pj)		
		for exp, act in zip(exp_ds.GetGeoTransform(), ds.GetGeoTransform()):
			self.assertAlmostEqual(exp, act, places=4)


class DateParsingTests(unittest.TestCase):

	def test_parse_short_date_pre2000(self):
		dstr = "980416"
		self.assertEqual(date(1998, 4, 16), roipac.parse_date(dstr))

	def test_parse_short_date_post2000(self):
		dstr = "081006"
		self.assertEqual(date(2008, 10, 6), roipac.parse_date(dstr))

	def test_parse_date_range(self):
		dstr = "980416-081006"
		exp = (date(1998, 4, 16), date(2008, 10, 6))
		self.assertEqual(exp, roipac.parse_date(dstr))


class HeaderParsingTests(unittest.TestCase):
	'''Verifies ROIPAC headers are parsed correctly.'''
	
	# short format header tests 
	def test_parse_short_roipac_header(self):
		hdrs = roipac.parse_header(SHORT_HEADER_PATH)
		self.assertEqual(hdrs[ifc.PYRATE_NCOLS], 47)
		self.assertEqual(hdrs[ifc.PYRATE_NROWS], 72)
		self.assertAlmostEqual(hdrs[ifc.PYRATE_LONG], 150.910)
		self.assertEqual(hdrs[ifc.PYRATE_X_STEP], 0.000833333)
		self.assertEqual(hdrs[ifc.PYRATE_LAT], -34.170000000)
		self.assertEqual(hdrs[ifc.PYRATE_Y_STEP], -0.000833333)
		self.assertEqual(hdrs[ifc.PYRATE_WAVELENGTH_METRES], 0.0562356424)

	def test_parse_short_header_has_timespan(self):
		# Ensures TIME_SPAN_YEAR field is added during parsing
		hdrs = roipac.parse_header(SHORT_HEADER_PATH)
		self.assertTrue(hdrs.has_key(roipac.TIME_SPAN_YEAR))

		# check time span calc
		master = date(2006, 06, 19)
		slave = date(2006, 10, 02)
		diff = (slave - master).days / 365.25
		self.assertEqual(diff, hdrs[roipac.TIME_SPAN_YEAR])


	# long format header tests
	def test_parse_full_roipac_header(self):
		# Ensures "long style" original header can be parsed correctly
		hdrs = roipac.parse_header(FULL_HEADER_PATH)

		# check DATE/ DATE12 fields are parsed correctly
		date0 = date(2006, 6, 19) # from "DATE 060619" header
		date2 = date(2006, 8, 28) # from DATE12 060619-060828
		self.assertEqual(hdrs[ifc.PYRATE_DATE], date0)
		self.assertEqual(hdrs[ifc.PYRATE_DATE2], date2)

	def test_read_full_roipac_header2(self):
		# Tests header from cropped original dataset is parsed correctly
		hdrs = roipac.parse_header(FULL_HEADER_PATH)
		self.assertTrue(len(hdrs) is not None)

	def test_xylast(self):
		# Test the X_LAST and Y_LAST header elements are calculated
		hdrs = roipac.parse_header(FULL_HEADER_PATH)
		self.assertAlmostEqual(hdrs[roipac.X_LAST], 151.8519444445)
		self.assertAlmostEqual(hdrs[roipac.Y_LAST], -34.625)
