'''
Created on 17/09/2012
@author: bpd900
'''

import unittest, datetime
from glob import glob
from os.path import join
from numpy.testing import assert_array_almost_equal

import config
from shared import Ifg
from ifgconstants import X_FIRST, Y_FIRST, WIDTH, FILE_LENGTH, X_STEP, Y_STEP


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
		def str2date(s):
			segs = s[:4], s[4:6], s[6:] # year, month, day
			return datetime.date(*[int(s) for s in segs])

		raw_date = ['20060619', '20060828', '20061002', '20061106', '20061211',
							'20070115', '20070219', '20070326', '20070430', '20070604',
							'20070709', '20070813', '20070917']
		exp_date = [str2date(d) for d in raw_date]
		exp_repeat = [1,1,3,3,4,3,3,3,3,3,3,2,2]

		exp_span = [0, 0.1916, 0.2875, 0.3833, 0.4791, 0.5749, 0.6708, 0.7666,
							0.8624, 0.9582, 1.0541, 1.1499, 1.2457]

		# test against Hua's results
		base = "../../tests/sydney_test/obs"
		paths = join(base, "ifms_17")
		with open(paths) as f:
			ifgs = [Ifg(join(base, path)) for path in f.readlines()]
			epochs = config.get_epochs(ifgs)

		self.assertEqual(len(exp_date), len(epochs.date) )
		self.assertTrue((exp_date == epochs.date).all())
		self.assertTrue((exp_repeat == epochs.repeat).all())
		assert_array_almost_equal(exp_span, epochs.span, decimal=4)


class IfgPrepTests(unittest.TestCase):

	def setUp(self):
		search = "../../tests/sydney_test/obs/*.unw"
		files = glob(search)
		self.ifgs = [Ifg(f) for f in files]
		self.assertTrue(len(self.ifgs) > 0)

	def test_check_xy_steps(self):
		self.assertTrue(config._check_xy_extents(self.ifgs))

	def test_check_xy_steps_failure(self):
		# test failure from breaking X and Y steps
		self.ifgs[1].X_STEP *= 2
		self.assertEqual(config._check_xy_extents(self.ifgs), X_STEP)
		self.ifgs[1].X_STEP /= 2 # revert X step
		self.ifgs[1].Y_STEP /= 2
		self.assertEqual(config._check_xy_extents(self.ifgs), Y_STEP)

	def test_check_xy_extents(self):
		self.assertTrue(config._check_xy_extents(self.ifgs))

	def test_check_xy_extents_failure(self):
		# test failure due to breaking X and Y extents
		self.ifgs[2].X_FIRST += 10
		self.assertEqual(config._check_xy_extents(self.ifgs), X_FIRST)
		self.ifgs[2].X_FIRST -= 10 # revert
		self.ifgs[2].Y_FIRST /= 2
		self.assertEqual(config._check_xy_extents(self.ifgs), Y_FIRST)

	def test_check_xy_extents_failure2(self):
		# validate length and width
		self.ifgs[3].WIDTH *= 2
		self.assertEqual(config._check_xy_extents(self.ifgs), WIDTH)
		self.ifgs[3].WIDTH /= 2 # revert X extent
		self.ifgs[3].FILE_LENGTH /= 2
		self.assertEqual(config._check_xy_extents(self.ifgs), FILE_LENGTH)

	def test_prepare_ifgs_no_looks(self):
		params = { config.IFG_LKSX : 0, config.IFG_LKSY : 0, config.IFG_CROP_OPT : None }
		config.prepare_ifgs(self.ifgs, params)

	def test_prepare_ifgs_looks(self):
		params = { config.IFG_LKSX : 3, config.IFG_LKSY : 3 }
		config.prepare_ifgs(self.ifgs, params)




if __name__ == "__main__":
	unittest.main()

