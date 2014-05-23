'''
Tests for geodesy/pyproj algorithms for PyRate.
Author:  Ben Davies
Created: 13/3/13
'''

import unittest

from pyrate.geodesy import cell_size, utm_zone


class GeodesyTests(unittest.TestCase):

	def test_utm_zone(self):
		# test some different zones
		for lon in [174.0, 176.5, 179.999, 180.0]:
			self.assertEqual(60, utm_zone(lon))

		for lon in [144.0, 144.1, 146.3456, 149.9999]:
			self.assertEqual(55, utm_zone(lon))

		for lon in [-180.0, -179.275, -176.925]:
			self.assertEqual(1, utm_zone(lon))

		for lon in [-72.0, -66.1]:
			self.assertEqual(19, utm_zone(lon))

		for lon in [0.0, 0.275, 3.925, 5.999]:
			self.assertEqual(31, utm_zone(lon))


	def test_cell_size_polar_region(self):
		# Can't have polar area zones: http://www.dmap.co.uk/utmworld.htm
		for lat in [-80.1, -85.0, -90.0, 84.1, 85.0, 89.9999, 90.0]:
			self.assertRaises(ValueError, cell_size, lat, 0, 0.1, 0.1)


	def test_cell_size_calc(self):
		# test conversion of X|Y_STEP to X|Y_SIZE
		x_deg = 0.000833333
		y_deg = -x_deg

		approx = 90.0 # x_deg is approx 90m
		exp_low = approx - (.15 * approx) # assumed tolerance
		exp_high = approx + (.15 * approx)

		latlons = [(10.0, 15.0), (-10.0, 15.0), (10.0, -15.0), (-10.0, -15.0),
		           (178.0, 33.0), (-178.0, 33.0), (178.0, -33.0), (-178.0, -33.0) ]
		for lon, lat in latlons:
			xs, ys = cell_size(lat, lon, x_deg, y_deg)
			for s in (xs, ys):
				self.assertTrue(s > 0, msg="size=%s" % s)
				self.assertTrue(s > exp_low, msg="size=%s" % s)
				self.assertTrue(s < exp_high, msg="size=%s" % s)
