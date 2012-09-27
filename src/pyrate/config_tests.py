'''
Created on 17/09/2012
@author: bpd900
'''

import unittest
import glob

from shared import Ifg
import config


class ConfigTests(unittest.TestCase):
	
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


	def test_get_epochs(self):
		exp_date = [20060619, 20060828, 20061002, 20061106, 20061211, 20070115,
							20070219, 20070326, 20070430, 20070604, 20070709, 20070813,
							20070917]
		exp_repeat = [1,1,3,3,4,3,3,3,3,3,3,2,2]
		exp_span = [0, 0.1916, 0.2875, 0.3833, 0.4791, 0.5749, 0.6708, 0.7666,
							0.8624, 0.9582, 1.0541, 1.1499, 1.2457]
		
		paths = "../../tests/sydney_test/obs/geo*.unw"
		ifgs = [Ifg(path) for path in glob.glob(paths)]
		epochs = config.get_epochs(ifgs)
		self.assertEqual(exp_date, epochs.date)
		self.assertEqual(exp_repeat, epochs.repeat)
		self.assertEqual(exp_span, epochs.span)


if __name__ == "__main__":
	unittest.main()