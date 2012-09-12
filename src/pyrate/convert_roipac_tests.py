'''
Created on 12/09/2012
@author: Ben Davies, ANUSF
				 ben.davies@anu.edu.au
'''


import unittest
import convert_roipac
from shared import IfgConstants


class ConversionTests(unittest.TestCase):
	
	HEADER_PATH = "../../tests/sydney_test/obs/geo_060619-061002.unw.rsc"
	

	def test_read_roipac_header(self):
		f = open(self.HEADER_PATH)
		header_text = f.read()
		f.close()
		
		hdrs = convert_roipac._read_roipac_header(header_text)
		self.assertTrue(hdrs is not None)
		self.assertEqual(hdrs[IfgConstants.WIDTH], 47)
		self.assertEqual(hdrs[IfgConstants.FILE_LENGTH], 72)
		self.assertAlmostEqual(hdrs[IfgConstants.X_FIRST], 150.910)


	def test_read_roipac_header_file(self):
		hdrs = convert_roipac._read_roipac_header(self.HEADER_PATH)
		self.assertEqual(hdrs[IfgConstants.X_STEP], 0.000833333)
		self.assertEqual(hdrs[IfgConstants.Y_FIRST], -34.170000000)
		self.assertEqual(hdrs[IfgConstants.Y_STEP], -0.000833333)
		self.assertEqual(hdrs[IfgConstants.WAVELENGTH], 0.0562356424)
