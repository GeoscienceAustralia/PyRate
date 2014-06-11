'''
Tests for the GAMMA file format interface

Created on 29/05/2014
@author: Ben Davies NCI
'''

import os
import sys
import unittest
from os.path import join
from datetime import date

from pyrate import gamma
import pyrate.ifgconstants as ifc
from pyrate.tests.common import GAMMA_TEST_DIR

from numpy.testing import assert_array_almost_equal

import gdal
gdal.UseExceptions()


LIGHTSPEED = 3e8 # approx


class GammaCommandLineTests(unittest.TestCase):

	def setUp(self):
		self.base = join(os.environ['PYRATEPATH'], 'tests', 'gamma')
		self.hdr = join(self.base, 'dem16x20raw.dem.par')

	def tearDown(self):
		os.remove(self.exp_path) 

	def test_cmd_ifg(self):
		data = join(self.base, '16x20_20090713-20090817_VV_4rlks_utm.unw')
		sys.argv = ['gamma.py', '-d', '/tmp', self.hdr, data]
		self.exp_path = '/tmp/16x20_20090713-20090817_VV_4rlks_utm.tif'
		self.common_check()

	def test_cmd_dem(self):
		data = join(self.base, 'dem16x20raw.dem')
		sys.argv = ['gamma.py', '-d', '/tmp', self.hdr, data]
		self.exp_path = '/tmp/dem16x20raw.tif'
		self.common_check()

	def common_check(self):
		gamma.main()
		self.assertTrue(os.path.exists(self.exp_path))


class GammaToGeoTiffTests(unittest.TestCase):
	'Tests conversion of GAMMA rasters to custom PyRate GeoTIFF'

	@classmethod
	def setUpClass(cls):
		# create common combined header obj so the headers are only read once 
		# tricker: needs both ifg headers, and DEM one for the extents 
		filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
		hdr_paths = [join(GAMMA_TEST_DIR, f) for f in filenames]
		hdrs = [gamma.parse_epoch_header(p) for p in hdr_paths]
		dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')

		cls.DEM_HDR = gamma.parse_dem_header(dem_hdr_path) 
		cls.COMBINED = gamma.combine_headers(*hdrs, dem_hdr=cls.DEM_HDR)

	def test_to_geotiff_dem(self):
		hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
		hdr = gamma.parse_dem_header(hdr_path)
		data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
		dest = "/tmp/tmp_gamma_dem.tif"

		gamma.to_geotiff(hdr, data_path, dest, nodata=0)
		exp_path = join(GAMMA_TEST_DIR, 'dem16x20_subset_from_gamma.tif')
		exp_ds = gdal.Open(exp_path)
		ds = gdal.Open(dest)

		# compare data and geographic headers
		assert_array_almost_equal(exp_ds.ReadAsArray(), ds.ReadAsArray())
		self.compare_rasters(ds, exp_ds)
		md = ds.GetMetadata()
		self.assertTrue(md['AREA_OR_POINT'] == 'Area')

	def test_to_geotiff_ifg(self):
		dest = '/tmp/tmp_gamma_ifg.tif'
		data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.unw')
		gamma.to_geotiff(self.COMBINED, data_path, dest, nodata=0)

		ds = gdal.Open(dest)
		exp_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.tif')
		exp_ds = gdal.Open(exp_path)

		# compare data and geographic headers
		assert_array_almost_equal(exp_ds.ReadAsArray(), ds.ReadAsArray())
		self.compare_rasters(ds, exp_ds)

		md = ds.GetMetadata()
		self.assertEqual(len(md), 5)
		self.assertTrue(md[ifc.PYRATE_DATE] == str(date(2009, 7, 13)))
		self.assertTrue(md[ifc.PYRATE_DATE2] == str(date(2009, 8, 17)))
		self.assertTrue(md[ifc.PYRATE_TIME_SPAN] == str(35 / ifc.DAYS_PER_YEAR))

		wavelen = float(md[ifc.PYRATE_WAVELENGTH_METRES])
		self.assertAlmostEqual(wavelen, 0.05627457792190739)

	def test_to_geotiff_wrong_input_data(self):
		# use TIF, not UNW for data
		dest = '/tmp/tmp_gamma_ifg.tif'
		data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.tif')
		self.assertRaises(gamma.GammaException, gamma.to_geotiff,
							self.COMBINED, data_path, dest, nodata=0)

	def test_mismatching_cell_resolution(self):
		hdrs = self.DEM_HDR.copy()
		hdrs[ifc.PYRATE_X_STEP] = 0.1 # fake a mismatch
		data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.unw')
		dest = '/tmp/fake'

		self.assertRaises(gamma.GammaException, gamma.to_geotiff, hdrs, data_path, dest, 0)

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

	def test_bad_projection(self):
		hdr = self.DEM_HDR.copy() 
		hdr[ifc.PYRATE_DATUM] = 'nonexistent projection'
		data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
		dest = "/tmp/tmp_gamma_dem2.tif"
		self.assertRaises(gamma.GammaException, gamma.to_geotiff, hdr, data_path, dest, nodata=0)


class GammaHeaderParsingTests(unittest.TestCase):
	'Tests conversion of GAMMA headers to Py dicts'

	def test_parse_gamma_epoch_header(self):
		# minimal required headers are:
		# date:      2009  7 13
		# radar_frequency:        5.3310040e+09   Hz
		path = join(GAMMA_TEST_DIR, 'r20090713_VV.slc.par')
		hdrs = gamma.parse_epoch_header(path)

		exp_date = date(2009, 7, 13)
		self.assertEqual(hdrs[ifc.PYRATE_DATE], exp_date)

		exp_wavelen = LIGHTSPEED / 5.3310040e+09
		self.assertEqual(hdrs[ifc.PYRATE_WAVELENGTH_METRES], exp_wavelen)


	def test_parse_gamma_dem_header(self):
		path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
		hdrs = gamma.parse_dem_header(path)

		self.assertEqual(hdrs[ifc.PYRATE_NCOLS], 16)
		self.assertEqual(hdrs[ifc.PYRATE_NROWS], 20)
		self.assertEqual(hdrs[ifc.PYRATE_LAT], -33.3831945)
		self.assertEqual(hdrs[ifc.PYRATE_LONG], 150.3870833)
		self.assertEqual(hdrs[ifc.PYRATE_X_STEP], 6.9444445e-05)
		self.assertEqual(hdrs[ifc.PYRATE_Y_STEP], -6.9444445e-05)


# Test data for the epoch header combination
H0 = { ifc.PYRATE_DATE : date(2009, 7, 13),
		ifc.PYRATE_WAVELENGTH_METRES : 1.8,
	}

H1 = { ifc.PYRATE_DATE : date(2009, 8, 17),
		ifc.PYRATE_WAVELENGTH_METRES : 1.8,
	}

H1_ERR = { ifc.PYRATE_DATE : date(2009, 8, 17),
			ifc.PYRATE_WAVELENGTH_METRES : 2.4,
	}


class HeaderCombinationTests(unittest.TestCase):
	'Tests GAMMA epoch and DEM headers can be combined into a single Py dict'

	def setUp(self):
		self.err = gamma.GammaException
		dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
		self.dh = gamma.parse_dem_header(dem_hdr_path)

	def test_combine_headers(self):
		filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
		paths = [join(GAMMA_TEST_DIR, p) for p in filenames]
		hdr0, hdr1 = [gamma.parse_epoch_header(p) for p in paths]

		chdr = gamma.combine_headers(hdr0, hdr1, self.dh)

		exp_timespan = (18 + 17) / ifc.DAYS_PER_YEAR
		self.assertEqual(chdr[ifc.PYRATE_TIME_SPAN], exp_timespan)

		exp_date = date(2009, 7, 13)
		self.assertEqual(chdr[ifc.PYRATE_DATE], exp_date)
		exp_date2 = date(2009, 8, 17)
		self.assertEqual(chdr[ifc.PYRATE_DATE2], exp_date2)

		exp_wavelen = LIGHTSPEED / 5.3310040e+09
		self.assertEqual(chdr[ifc.PYRATE_WAVELENGTH_METRES], exp_wavelen)

	def test_fail_non_dict_header(self):
		self.assertRaises(self.err, gamma.combine_headers, H0, '', self.dh)
		self.assertRaises(self.err, gamma.combine_headers, '', H0, self.dh)
		self.assertRaises(self.err, gamma.combine_headers, H0, H1, None)
		self.assertRaises(self.err, gamma.combine_headers, H0, H1, '')

	def test_fail_mismatching_wavelength(self):
		self.assertRaises(self.err, gamma.combine_headers, H0, H1_ERR, self.dh)

	def test_fail_same_date(self):
		self.assertRaises(self.err, gamma.combine_headers, H0, H0, self.dh)

	def test_fail_bad_date_order(self):
		self.assertRaises(self.err, gamma.combine_headers, H1, H0, self.dh)
