'''
Tests for the GAMMA file format interface

Created on 29/05/2014
@author: Ben Davies NCI
'''

import unittest
from os.path import join

from datetime import date
from pyrate import gamma
from pyrate.tests.common import HEADERS_TEST_DIR, GAMMA_TEST_DIR

from numpy.testing import assert_array_almost_equal

import gdal
gdal.UseExceptions()


LIGHTSPEED = 3e8 # approx



# TODO: needs a chopped out segment of data (ie, 16x20 pixels)
class GammaToGeoTiffTests(unittest.TestCase):
	
	def test_to_geotiff_dem(self):
		# TODO: refactor with Ifg class later
		hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
		hdr = gamma.parse_dem_header(hdr_path)
		data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
		dest = "/tmp/tmpdem.tif"
		
		# TODO: refactor to take a header path (or autodetect header)
		gamma.to_geotiff(hdr, data_path, dest, nodata=0)
		ds = gdal.Open(dest)
		band = ds.GetRasterBand(1)
		data = band.ReadAsArray()
		
		exp_path = join(GAMMA_TEST_DIR, 'dem16x20_subset_from_gamma.tif')
		exp_ds = gdal.Open(exp_path)
		exp_band = exp_ds.GetRasterBand(1)
		exp_data = exp_band.ReadAsArray() 

		# compare data and metadata
		assert_array_almost_equal(exp_data, data)
		nodata = band.GetNoDataValue()
		self.assertFalse(nodata is None)
		self.assertEqual(exp_band.GetNoDataValue(), nodata)
		
		self.assertEqual(exp_ds.GetProjection(), ds.GetProjection())		
		for exp, act in zip(exp_ds.GetGeoTransform(), ds.GetGeoTransform()):
			self.assertAlmostEqual(exp, act, places=4)
	
	def test_to_geotiff_ifg(self):
		raise NotImplementedError


class GammaHeaderParsingTests(unittest.TestCase):
	'TODO'

	def test_parse_gamma_epoch_header(self):
		# minimal required headers are:
		# date:      2009  7 13
		# radar_frequency:        5.3310040e+09   Hz
		path = join(HEADERS_TEST_DIR, 'r20090713_VV.slc.par')
		hdrs = gamma.parse_epoch_header(path)
		
		exp_date = date(2009, 7, 13)
		self.assertEqual(hdrs['DATE'], exp_date)
		
		exp_wavelen = LIGHTSPEED / 5.3310040e+09
		self.assertEqual(hdrs['WAVELENGTH_METRES'], exp_wavelen)


	def test_parse_gamma_dem_header(self):
		path = join(HEADERS_TEST_DIR, '20090713_VV_4rlks_utm_dem.par')
		hdrs = gamma.parse_dem_header(path)
		
		self.assertEqual(hdrs['NCOLS'], 20316)
		self.assertEqual(hdrs['NROWS'], 16872)
		self.assertEqual(hdrs['LAT'], -33.3831945)
		self.assertEqual(hdrs['LONG'], 150.3870833)
		self.assertEqual(hdrs['X_STEP'], 6.9444445e-05)
		self.assertEqual(hdrs['Y_STEP'], -6.9444445e-05)
	

# Test data for the epoch header combination
H0 = { 'DATE' : date(2009, 7, 13),
		'WAVELENGTH_METRES' : 1.8,
	}

H1 = { 'DATE' : date(2009, 8, 17),
		'WAVELENGTH_METRES' : 1.8,
	}

H1_FAULT = { 'DATE' : date(2009, 8, 17),
			'WAVELENGTH_METRES' : 2.4,
	}


class HeaderCombinationTests(unittest.TestCase):
	
	def test_combine_headers(self):
		filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
		paths = [join(HEADERS_TEST_DIR, p) for p in filenames]
		hdr0, hdr1 = [gamma.parse_epoch_header(p) for p in paths]
		chdr = gamma.combine_headers(hdr0, hdr1)
		
		exp_timespan = (18 + 17) / 365.25 
		self.assertEqual(chdr['TIME_SPAN_YEAR'], exp_timespan)
		
		# TODO: what to do with second date? Ignore as per ROIPAC?
		exp_date = date(2009, 7, 13)
		self.assertEqual(chdr['DATE'], exp_date)
		
		exp_wavelen = LIGHTSPEED / 5.3310040e+09
		self.assertEqual(chdr['WAVELENGTH_METRES'], exp_wavelen)
	
	def test_fail_mismatching_wavelength(self):
		self.assertRaises(gamma.GammaError, gamma.combine_headers, H0, H1_FAULT)
	
	def test_fail_same_date(self):
		self.assertRaises(gamma.GammaError, gamma.combine_headers, H0, H0)
		
	def test_fail_bad_date_order(self):
		self.assertRaises(gamma.GammaError, gamma.combine_headers, H1, H0)
