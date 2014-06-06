'''
Tests for the pyrate configuration file parser.

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import sys, unittest
from os.path import exists, join

from pyrate import config
from common import SYD_TEST_DIR, SYD_TEST_TIF


class ConfigTests(unittest.TestCase):

	def __init__(self, *args, **kwargs):
		super(ConfigTests, self).__init__(*args, **kwargs)
		if not exists(SYD_TEST_TIF):
			sys.exit("ERROR: Missing sydney_test data for unit tests\n")

	# TODO: are more conf parsing tests needed?

	def test_read_param_file(self):
		conf_file = join(SYD_TEST_DIR, 'pyrate.conf')
		params = config.parse_conf_file(conf_file)
		for k in params.keys():
			self.assertTrue(k != "")
			self.assertTrue(params[k] != "")
			self.assertFalse(k.endswith(":")) # are the colons removed?
	
	
	def test_read_param_file_no_refpixel(self):
		# ensure the parser can handle empty fields
		conf_file = join(SYD_TEST_DIR, 'pyrate2.conf')
		params = config.parse_conf_file(conf_file)
		
		self.assertTrue(params[config.REFX] == -1)
		self.assertTrue(params[config.REFY] == -1)


	def test_parse_namelist(self):
		nl = join(SYD_TEST_TIF, 'ifms_17')
		result = config.parse_namelist(nl)
		self.assertEqual(len(result), 17)
		files = ["geo_060619-061002.tif", "geo_060828-061211.tif",
					"geo_061002-070430.tif", "geo_070115-070917.tif",
					"geo_070219-070604.tif" ]

		for path in files:
			self.assertTrue(path in result)


if __name__ == "__main__":
	unittest.main()
