'''
Created on 12/09/2012
@author: Ben Davies, ANUSF
				 ben.davies@anu.edu.au
'''


import os, unittest
import datetime
from osgeo import gdal

import roipac
import ifgconstants as IFC


class ConversionTests(unittest.TestCase):
	
	SHORT_HEADER_PATH = "../../tests/sydney_test/obs/geo_060619-061002.unw.rsc"
	FULL_HEADER_PATH  = "../../tests/headers/geo_060619-060828.unw.rsc"
	FULL_HEADER_PATH2 = "../../tests/single/geo_060619-061002.unw.rsc"
	

	def test_read_short_roipac_header(self):
		f = open(self.SHORT_HEADER_PATH)
		header_text = f.read()
		f.close()
		
		hdrs = roipac.parse_header(header_text)
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


	def test_ehdr_header_defaults(self):
		# test default header filename
		base_hdr = os.path.abspath("../../tests/single/geo_060619-061002.unw.rsc")
		hdr = "/tmp/geo_060619-061002.unw.rsc"
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
		gdal.UseExceptions()
		ds = gdal.Open(exp_data)
		self.assertTrue(ds is not None)
		bands = ds.GetRasterBand(1), ds.GetRasterBand(1)
		self.assertTrue(all(bands))
		
		del bands, ds
		os.unlink(exp_data)
		os.remove(exp_hdr)
		os.unlink(hdr)


	def test_convert_roipac(self):
		raise NotImplementedError
		#expected = "../../tests/sydney_test/obs/geo_060619-061002.tif"
		#roipac.roipac(src, dest, fmt)
