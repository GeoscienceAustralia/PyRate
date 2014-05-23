'''
Created on 17/09/2012
@author: bpd900
'''

import sys, unittest
from os.path import exists, join

import config
from tests_common import SYD_TEST_DIR, SYD_TEST_OBS


class ConfigTests(unittest.TestCase):

	def __init__(self, *args, **kwargs):
		super(ConfigTests, self).__init__(*args, **kwargs)
		if not exists(SYD_TEST_OBS):
			sys.exit("ERROR: Missing sydney_test data for unit tests\n")

	# TODO: are more conf parsing tests needed?

	def test_read_param_file(self):
		conf_file = join(SYD_TEST_DIR, 'pyrate.conf')
		params = config.parse_conf_file(conf_file)  # (params_file)
		for k in params.keys():
			self.assertTrue(k != "")
			self.assertTrue(params[k] != "")
			self.assertFalse(k.endswith(":")) # have we gotten rid of the colons?


	def test_parse_namelist(self):
		nl = join(SYD_TEST_OBS, 'ifms_17')
		result = config.parse_namelist(nl)
		self.assertEqual(len(result), 17)
		files = ["geo_060619-061002.unw", "geo_060828-061211.unw",
					"geo_061002-070430.unw", "geo_070115-070917.unw",
					"geo_070219-070604.unw" ]

		for path in files:
			self.assertTrue(path in result)


if __name__ == "__main__":
	unittest.main()
