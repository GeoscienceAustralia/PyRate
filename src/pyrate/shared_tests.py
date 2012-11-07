'''
Created on 12/09/2012
@author: bpd900
'''

import unittest, os
from itertools import product
from numpy import isnan

from gdal import Dataset, UseExceptions
UseExceptions()

import shared
from shared import Ifg, IfgException, Raster, RasterException
from ifgconstants import Z_OFFSET, Z_SCALE, PROJECTION, DATUM



class IfgTests(unittest.TestCase):
	'''Unit tests for the Ifg/interferogram class.'''

	def setUp(self):
		self.ifg = Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')


	def test_create_ifg(self):
		self.assertTrue(os.path.exists(self.ifg.hdr_path)) # validate header path


	def test_headers_as_attr(self):
		attrs = ['WIDTH', 'FILE_LENGTH', 'X_FIRST', 'X_STEP',
			'Y_FIRST', 'Y_STEP', 'WAVELENGTH', 'MASTER', 'SLAVE']

		for a in attrs:
			self.assertTrue(getattr(self.ifg, a) is not None)


	def test_open(self):
		self.assertTrue(self.ifg.dataset is None)
		self.ifg.open()
		self.assertTrue(self.ifg.dataset is not None)
		self.assertTrue(isinstance(self.ifg.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(IfgException, self.ifg.open)
		os.remove(self.ifg.ehdr_path)


	def test_xylast(self):
		# ensure the X|Y_LAST header element has been created
		self.ifg.open()
		self.assertAlmostEqual(self.ifg.X_LAST, 150.9491667)
		self.assertAlmostEqual(self.ifg.Y_LAST, -34.23)


	def test_num_cells(self):
		self.ifg.open()
		data = self.ifg.amp_band.ReadAsArray()
		ys, xs = data.shape
		exp_ncells = ys * xs
		self.assertEqual(exp_ncells, self.ifg.num_cells)


	def test_amp_band(self):
		try:
			_ = self.ifg.amp_band
			self.fail("Should not be able to access band without open dataset")
		except IfgException, ex:
			pass

		self.ifg.open()
		data = self.ifg.amp_band.ReadAsArray()
		self.assertEqual(data.shape, (72,47) )


	def test_phase_band(self):
		try:
			_ = self.ifg.phase_band
			self.fail("Should not be able to access band without open dataset")
		except IfgException, ex:
			pass

		self.ifg.open()
		data = self.ifg.phase_band.ReadAsArray()
		self.assertEqual(data.shape, (72,47) )


	def test_nan_fraction(self):
		try:
			# NB: self.assertRaises doesn't work here
			_ = self.ifg.nan_fraction()
			self.fail("Shouldn't be able to call nan_craction() with unopened Ifg")
		except IfgException as ex:
			pass

		self.ifg.open()
		data = self.ifg.dataset.GetRasterBand(2).ReadAsArray()
		ys, xs = data.shape

		nans = 0
		for y,x in product(xrange(ys), xrange(xs)):
			if isnan(data[y,x]): nans += 1

		num_cells = float(ys * xs)
		nan_frac = nans / num_cells
		self.assertEqual(nan_frac, self.ifg.nan_fraction)



class RasterTests(unittest.TestCase):
	'''Unit tests for the generic Raster class.'''

	def setUp(self):
		self.ras = Raster('../../tests/sydney_test/dem/sydney_trimmed.dem')


	def test_create_raster(self):
		self.assertTrue(os.path.exists(self.ras.hdr_path)) # validate header path


	def test_headers_as_attr(self):
		attrs = ['WIDTH', 'FILE_LENGTH', 'X_FIRST', 'X_STEP',
			'Y_FIRST', 'Y_STEP', Z_OFFSET, Z_SCALE, PROJECTION, DATUM]

		for a in attrs:
			self.assertTrue(getattr(self.ras, a) is not None)


	def test_is_dem(self):
		self.assertTrue(self.ras.is_dem())
		self.ras = Raster('../../tests/sydney_test/obs/geo_060619-061002.unw')
		self.assertFalse(self.ras.is_dem())


	def test_open(self):
		self.assertTrue(self.ras.dataset is None)
		self.ras.open()
		self.assertTrue(self.ras.dataset is not None)
		self.assertTrue(isinstance(self.ras.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(RasterException, self.ras.open)
		os.remove(self.ras.ehdr_path)


	def test_band(self):
		try:
			_ = self.ras.band
			self.fail("Should not be able to access band without open dataset")
		except RasterException, ex:
			pass

		self.ras.open()
		data = self.ras.band.ReadAsArray()
		self.assertEqual(data.shape, (72,47) )


if __name__ == "__main__":
	unittest.main()
