'''
Tests for the GAMMA file format interface

Created on 29/05/2014
@author: Ben Davies NCI
'''

import unittest
from os.path import join

from datetime import date
from pyrate import gamma
from pyrate.tests.common import HEADERS_TEST_DIR


# TODO: test class for funcs dealing combining 2 headers
# TODO: create a time_span_year item (needs 2 headers)

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
		
		speed_of_light = 3e8 # approx
		exp_wavelen = speed_of_light / 5.3310040e+09
		self.assertEqual(hdrs['WAVELENGTH_METRES'], exp_wavelen)


	def test_parse_gamma_dem_header(self):
		path = join(HEADERS_TEST_DIR, '20090713_VV_4rlks_utm_dem.par')
		hdrs = gamma.parse_dem_header(path)
		
		self.assertEqual(hdrs['NCOLS'], 20316)
		self.assertEqual(hdrs['NROWS'], 16872)
		self.assertEqual(hdrs['LAT'], -33.3831945)
		self.assertEqual(hdrs['LONG'], 150.3870833)
		self.assertEqual(hdrs['X_STEP'], 6.9444445e-05)
		self.assertEqual(hdrs['Y_STEP'], 6.9444445e-05)



class HeaderCombinationTests(unittest.TestCase):
	pass

	# TODO:
	# order the dates correctly
	# dates are not the same
	# frequencies/wavelens match
	# 
		
		
		
