'''
Created on 12/09/2012
@author: Ben Davies, ANUSF
         ben.davies@anu.edu.au
'''


import os, sys
from os.path import exists
import unittest, datetime

from numpy import amin, amax, zeros
from numpy.testing import assert_array_equal

from gdal import Open, UseExceptions
UseExceptions()

import roipac
from roipac import RoipacException
import ifgconstants as IFC


class ConversionTests(unittest.TestCase):
	'''Verifies conversion of ROIPAC files to EHdr format.'''

	def __init__(self, *args, **kwargs):
		super(ConversionTests, self).__init__(*args, **kwargs)
		if not exists("../../tests/headers"):
			sys.exit("ERROR: Missing the 'headers' data for unittests\n")

	SHORT_HEADER_PATH = "../../tests/sydney_test/obs/geo_060619-061002.unw.rsc"
	FULL_HEADER_PATH  = "../../tests/headers/geo_060619-060828.unw.rsc"
	FULL_HEADER_PATH2 = "../../tests/single/geo_060619-061002.unw.rsc"

	def test_read_short_roipac_header(self):
		hdrs = roipac.parse_header(self.SHORT_HEADER_PATH)
		self.assertTrue(hdrs is not None)
		self.assertEqual(hdrs[IFC.WIDTH], 47)
		self.assertEqual(hdrs[IFC.FILE_LENGTH], 72)
		self.assertAlmostEqual(hdrs[IFC.X_FIRST], 150.910)


	def test_read_roipac_header_file(self):
		hdrs = roipac.parse_header(self.SHORT_HEADER_PATH)
		self.assertEqual(hdrs[IFC.X_STEP], 0.000833333)
		self.assertEqual(hdrs[IFC.Y_FIRST], -34.170000000)
		self.assertEqual(hdrs[IFC.Y_STEP], -0.000833333)
		self.assertEqual(hdrs[IFC.WAVELENGTH], 0.0562356424)


	def test_parse_short_date_pre2000(self):
		dstr = "980416"
		self.assertEqual(datetime.date(1998, 4, 16), roipac.parse_date(dstr))


	def test_parse_short_date_post2000(self):
		dstr = "081006"
		self.assertEqual(datetime.date(2008, 10, 6), roipac.parse_date(dstr))


	def test_parse_date_range(self):
		dstr = "980416-081006"
		exp = (datetime.date(1998, 4, 16), datetime.date(2008, 10, 6))
		self.assertEqual(exp, roipac.parse_date(dstr))


	def test_read_full_roipac_header(self):
		"Tests full original header can be parsed correctly"
		hdrs = roipac.parse_header(self.FULL_HEADER_PATH)

		# check dates are parsed correctly
		date0 = datetime.date(2006, 6, 19) # from  "DATE 060619" header
		date12 = (date0, datetime.date(2006, 8, 28)) # from DATE12   060619-060828
		self.assertEqual(hdrs[IFC.DATE], date0)
		self.assertEqual(hdrs[IFC.DATE12], date12)


	def test_date_alias(self):
		"""Test header has MASTER and SLAVE dates as keys"""
		hdrs = roipac.parse_header(self.FULL_HEADER_PATH)
		self.assertTrue(hdrs.has_key(IFC.MASTER))
		self.assertTrue(hdrs.has_key(IFC.SLAVE))
		self.assertEqual(hdrs[IFC.DATE], hdrs[IFC.MASTER])
		self.assertEqual(hdrs[IFC.DATE12][-1], hdrs[IFC.SLAVE])


	def test_timespan(self):
		"""Ensures the TIME_SPAN_YEAR element is present after parsing short header"""
		hdrs = roipac.parse_header(self.SHORT_HEADER_PATH)
		self.assertTrue(hdrs.has_key(IFC.TIME_SPAN_YEAR))

		# check time span calc
		master = datetime.date(2006, 06, 19)
		slave = datetime.date(2006, 10, 02)
		diff = (slave - master).days / 365.25
		self.assertEqual(diff, hdrs[IFC.TIME_SPAN_YEAR])


	def test_xylast(self):
		# Test the X_LAST and Y_LAST header elements are added
		hdrs = roipac.parse_header(self.FULL_HEADER_PATH)
		self.assertAlmostEqual(hdrs[IFC.X_LAST], 151.8519444445)
		self.assertAlmostEqual(hdrs[IFC.Y_LAST], -34.625)


	def test_read_full_roipac_header2(self):
		"Tests header from cropped original dataset is parsed correctly"
		hdrs = roipac.parse_header(self.FULL_HEADER_PATH)
		self.assertTrue(hdrs is not None)


	def test_filename_pair(self):
		base = "project/run/geo_070709-080310.unw"
		exp_hdr = "project/run/geo_070709-080310.unw.rsc"
		result = roipac.filename_pair(base)
		self.assertEqual( (base, exp_hdr), result)


	def test_roipac_to_ehdr_header(self):
		dest = "/tmp/ehdr.hdr"
		hdr = "../../tests/headers/geo_060619-060828.unw.rsc"
		roipac.to_ehdr_header(hdr, dest)

		with open(dest) as f:
			text = f.read()
			self.assertTrue("ncols 5451" in text)
			self.assertTrue("nrows 4366" in text)
			self.assertTrue("cellsize 0.0002777" in text) # compare to 7 places
			self.assertTrue("xllcorner 150.3377777" in text) # compare to 7 places
			exp_yll = "yllcorner " + str(round(-33.4122222 - (4366 * 0.0002777), 3))
			self.assertTrue(exp_yll in text, "Got " + exp_yll)

		os.remove(dest)


	def test_to_ehdr_header_defaults(self):
		# test default header filename
		base_hdr = os.path.abspath("../../tests/single/geo_060619-061002.unw.rsc")
		hdr = "/tmp/geo_060619-061002.unw.rsc"
		if os.path.exists(hdr):
			os.unlink(hdr)
		os.symlink(base_hdr, hdr)
		exp_hdr = "/tmp/geo_060619-061002.hdr"

		if os.path.exists(exp_hdr): os.remove(exp_hdr)
		roipac.to_ehdr_header(hdr)
		self.assertTrue(os.path.exists(exp_hdr))

		# add data to /tmp for GDAL test
		base_data = os.path.abspath("../../tests/single/geo_060619-061002.unw")
		exp_data = "/tmp/geo_060619-061002.unw"
		os.symlink(base_data, exp_data)

		# test GDAL can open the data with the new header & cleanup
		ds = Open(exp_data)
		self.assertTrue(ds is not None)
		bands = ds.GetRasterBand(1), ds.GetRasterBand(1)
		self.assertTrue(all(bands)) # both bands exist?

		del bands, ds
		os.unlink(exp_data)
		os.remove(exp_hdr)
		os.unlink(hdr)


	def test_to_ehdr_header_with_data(self):
		# ensure giving the data file breaks to_ehdr_header()
		src = "../../tests/sydney_test/obs/geo_060619-061002.unw"
		try:
			roipac.to_ehdr_header(src)
			self.fail("Should not be able to accept .unw data file")
		except:
			pass


	def test_to_ehdr_header_with_missing_file(self):
		# ensure giving the data file breaks to_ehdr_header()
		src = "../../tests/sydney_test/obs/fake.unw.rsc"
		self.assertRaises(IOError, roipac.to_ehdr_header, src)


	def test_to_ehdr_header_with_dir(self):
		# ensure giving the data file breaks to_ehdr_header()
		src = "../../tests/sydney_test/obs"
		self.assertRaises(IOError, roipac.to_ehdr_header, src)


	def test_to_ehdr_header_with_dem(self):
		dem_hdr = "../../tests/sydney_test/dem/sydney_trimmed.dem.rsc"
		act = roipac.to_ehdr_header(dem_hdr)
		self.assertEqual(act, dem_hdr[:-7] + "hdr")

		with open(act) as f:
			lines = [line.strip() for line in f.readlines()]
			values = [line.split() for line in lines]

		self.assertTrue(['ncols', '47'] in values)
		self.assertTrue(['nrows', '72'] in values)

		# TODO: test different forms of cellsize
		self.assertTrue(['cellsize', '0.000833333'] in values)
		#self.assertTrue(['xdim', '0.000277777777782'] in values)
		#self.assertTrue(['ydim', '-0.000277777777782'] in values)

		self.assertFalse(['nodata', '0'] in values)
		self.assertFalse(['nbands', '1'] in values)
		self.assertTrue(['byteorder', 'lsb'] in values)
		self.assertFalse(['layout', 'bil'] in values)
		self.assertTrue(['nbits', '16'] in values)
		self.assertTrue(['pixeltype', 'signedint'] in values)
		os.remove(act)


	def test_to_ehdr_header_gdal(self):
		# test data files can be opened, new headers generated, data is readable
		hdr = "../../tests/sydney_test/obs/geo_060619-061002.unw.rsc"
		ehdr = "../../tests/sydney_test/obs/geo_060619-061002.hdr"
		if os.path.exists(ehdr):
			os.remove(ehdr) # can be left behind if the test fails

		roipac.to_ehdr_header(hdr)
		self.assertTrue(os.path.exists(ehdr))

		# open with GDAL and ensure there is data
		src = "../../tests/sydney_test/obs/geo_060619-061002.unw"
		ds = Open(src)
		self.assertTrue(ds is not None)

		# check faked amplitude band 1 for 0s
		band = ds.GetRasterBand(1)
		data = band.ReadAsArray()
		shape = data.shape
		assert_array_equal(data, zeros(shape))

		# check phase (band 2) for data
		band = ds.GetRasterBand(2)
		nodata = band.GetNoDataValue()
		self.assertEqual(nodata, 0) # check default ROIPAC NODATA

		# ensure decent data is retrieved
		data = band.ReadAsArray()
		self.assertTrue(amin(data) != 0) # ignore max as NODATA is 0
		self.assertTrue(data.ptp() != 0)

		# cleanup
		os.remove(ehdr)
