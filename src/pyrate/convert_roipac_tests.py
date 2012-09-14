'''
Created on 12/09/2012
@author: Ben Davies, ANUSF
				 ben.davies@anu.edu.au
'''


import os, unittest
import convert_roipac
import ifgconstants as IFC


class ConversionTests(unittest.TestCase):
	
	HEADER_PATH = "../../tests/sydney_test/obs/geo_060619-061002.unw.rsc"
	

	def test_read_roipac_header(self):
		f = open(self.HEADER_PATH)
		header_text = f.read()
		f.close()
		
		hdrs = convert_roipac._read_roipac_header(header_text)
		self.assertTrue(hdrs is not None)
		self.assertEqual(hdrs[IFC.WIDTH], 47)
		self.assertEqual(hdrs[IFC.FILE_LENGTH], 72)
		self.assertAlmostEqual(hdrs[IFC.X_FIRST], 150.910)


	def test_read_roipac_header_file(self):
		hdrs = convert_roipac._read_roipac_header(self.HEADER_PATH)
		self.assertEqual(hdrs[IFC.X_STEP], 0.000833333)
		self.assertEqual(hdrs[IFC.Y_FIRST], -34.170000000)
		self.assertEqual(hdrs[IFC.Y_STEP], -0.000833333)
		self.assertEqual(hdrs[IFC.WAVELENGTH], 0.0562356424)
		
	
	def test_filename_pair(self):
		base = "project/run/geo_070709-080310.unw"
		exp_hdr = "project/run/geo_070709-080310.unw.rsc"
		result = convert_roipac.filename_pair(base)
		self.assertEqual( (base, exp_hdr), result)

		
	def test_roipac_to_ehdr_header(self):
		dest = "/tmp/ehdr.hdr"
		hdr = "../../tests/headers/geo_060619-060828.unw.rsc"
		convert_roipac.roipac_to_ehdr_header(hdr, dest)
		
		with open(dest) as f:
			text = f.read()
			self.assertTrue("ncols 5451" in text)
			self.assertTrue("nrows 4366" in text)
			self.assertTrue("cellsize 0.0002777" in text) # compare to 7 places
			self.assertTrue("xllcorner 150.3377777" in text) # compare to 7 places
			exp_yll = "yllcorner " + str(round(-33.4122222 - (4366 * 0.0002777), 3))
			self.assertTrue(exp_yll in text, "Got " + exp_yll)
			
		os.remove(dest)


	def test_roipac_to_ehdr_header_gdal_compatibility(self):
		# TODO: test default header filename is recognised by GDAL 
		raise NotImplementedError



	def test_convert_roipac(self):
		raise NotImplementedError
		#expected = "../../tests/sydney_test/obs/geo_060619-061002.tif"
		#convert_roipac.convert_roipac(src, dest, fmt)
		
	# TODO: def test_convert_roipac_dual_band(self):
	
