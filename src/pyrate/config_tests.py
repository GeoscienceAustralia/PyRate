'''
Created on 17/09/2012
@author: bpd900
'''

import sys, unittest
from os.path import exists

import config



class ConfigTests(unittest.TestCase):

	def __init__(self, *args, **kwargs):
		super(ConfigTests, self).__init__(*args, **kwargs)
		if not exists("../../tests/sydney_test/obs"):
			sys.exit("ERROR: Missing the sydney_test data for unittests\n")

	# TODO: are more conf parsing tests needed?

	def test_read_param_file(self):
		conf_file = "../../tests/sydney_test/pyrate.conf"
		params = config.parse_conf_file(conf_file)  # (params_file)
		for k in params.keys():
			self.assertTrue(k != "")
			self.assertTrue(params[k] != "")
			self.assertFalse(k.endswith(":")) # have we gotten rid of the colons?


	def test_parse_namelist(self):
		nl = "../../tests/sydney_test/obs/ifms_17"
		result = config.parse_namelist(nl)
		self.assertEqual(len(result), 17)
		files = ["geo_060619-061002.unw", "geo_060828-061211.unw",
						"geo_061002-070430.unw", "geo_070115-070917.unw",
						"geo_070219-070604.unw" ]

		for path in files:
			self.assertTrue(path in result)


if __name__ == "__main__":
	unittest.main()
