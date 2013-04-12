'''
Created on 12/09/2012
@author: bpd900
'''

import unittest, os
from itertools import product
from numpy import isnan, where, nan
from numpy.testing import assert_array_equal

from gdal import Open, Dataset, UseExceptions
UseExceptions()

from shared import Ifg, DEM, RasterException
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
		self.assertTrue(self.ifg.is_open is False)
		self.ifg.open()
		self.assertTrue(self.ifg.dataset is not None)
		self.assertTrue(self.ifg.is_open is True)
		self.assertTrue(isinstance(self.ifg.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(RasterException, self.ifg.open)
		os.remove(self.ifg.ehdr_path)


	def test_xylast(self):
		# ensure the X|Y_LAST header element has been created
		self.ifg.open()
		self.assertAlmostEqual(self.ifg.X_LAST, 150.9491667)
		self.assertAlmostEqual(self.ifg.Y_LAST, -34.23)


	def test_num_cells(self):
		# test cell size from header elements
		self.assertEqual(self.ifg.num_cells, 47 * 72)

		self.ifg.open()
		data = self.ifg.amp_band.ReadAsArray()
		ys, xs = data.shape
		exp_ncells = ys * xs
		self.assertEqual(exp_ncells, self.ifg.num_cells)


	def test_amp_band(self):
		try:
			_ = self.ifg.amp_band
			self.fail("Should not be able to access band without open dataset")
		except RasterException:
			pass

		self.ifg.open()
		data = self.ifg.amp_band.ReadAsArray()
		self.assertEqual(data.shape, (72,47) )


	def test_phase_band(self):
		try:
			_ = self.ifg.phase_band
			self.fail("Should not be able to access band without open dataset")
		except RasterException:
			pass

		self.ifg.open()
		data = self.ifg.phase_band.ReadAsArray()
		self.assertEqual(data.shape, (72,47) )


	def test_nan_count(self):
		self.ifg.open()
		num_nan = 0
		for row in self.ifg.phase_data:
			for v in row:
				if isnan(v):
					num_nan += 1

		self.assertEqual(num_nan, self.ifg.nan_count)


	def test_nan_fraction(self):
		try:
			# NB: self.assertRaises doesn't work here (as it is a property?)
			_ = self.ifg.nan_fraction
			self.fail("Shouldn't be able to call nan_fraction() with unopened Ifg")
		except RasterException:
			pass

		# NB: source data lacks 0 -> NaN conversion
		self.ifg.open()
		data = self.ifg.dataset.GetRasterBand(2).ReadAsArray()
		data = where(data == 0, nan, data) # fake 0 -> nan for the count below

		# manually count # nan cells
		nans = 0
		ys, xs = data.shape
		for y,x in product(xrange(ys), xrange(xs)):
			if isnan(data[y,x]): nans += 1
		del data

		num_cells = float(ys * xs)
		self.assertTrue(nans > 0)
		self.assertTrue(nans <= num_cells)
		self.assertEqual(nans / num_cells, self.ifg.nan_fraction)


	def test_phase_data_properties(self):
		# Use raw GDAL to isolate raster reading from Ifg functionality
		ds = Open(self.ifg.data_path)
		data = ds.GetRasterBand(2).ReadAsArray()
		del ds

		self.ifg.open()

		# test full array and row by row access
		assert_array_equal(data, self.ifg.phase_data)
		for y, row in enumerate(self.ifg.phase_rows):
			assert_array_equal(data[y], row)

		# test the data is cached if changed
		crd = (5,4)
		orig = self.ifg.phase_data[crd]
		self.ifg.phase_data[crd] *= 2
		nv = self.ifg.phase_data[crd] # pull new value out again
		self.assertEqual(nv, 2 * orig)


	def test_xy_size(self):
		self.ifg.open()
		self.assertFalse(self.ifg.X_SIZE is None)
		self.assertFalse(self.ifg.Y_SIZE is None)

		# test with tolerance from base 90m cell
		self.assertTrue(self.ifg.Y_SIZE > 88.0) # within 2% of cells over sydney?
		self.assertTrue(self.ifg.Y_SIZE < 92.0)

		syd_width = 76.9 # from nearby pirate coords
		self.assertTrue(self.ifg.X_SIZE > 0.97 * syd_width) # ~3% tolerance
		self.assertTrue(self.ifg.X_SIZE < 1.03 * syd_width)


	def test_centre_latlong(self):
		ifg = self.ifg
		ifg.open()
		lat_exp = ifg.Y_FIRST + ((ifg.FILE_LENGTH/2) * ifg.Y_STEP)
		long_exp = ifg.X_FIRST + ((ifg.WIDTH/2) * ifg.X_STEP)
		self.assertEqual(lat_exp, self.ifg.LAT_CENTRE)
		self.assertEqual(long_exp, self.ifg.LONG_CENTRE)


	def test_centre_cell(self):
		self.assertEqual(self.ifg.X_CENTRE, 23)
		self.assertEqual(self.ifg.Y_CENTRE, 36)



class DEMTests(unittest.TestCase):
	'''Unit tests for the generic Raster class.'''

	def setUp(self):
		self.ras = DEM('../../tests/sydney_test/dem/sydney_trimmed.dem')


	def test_create_raster(self):
		self.assertTrue(os.path.exists(self.ras.hdr_path)) # validate header path


	def test_headers_as_attr(self):
		attrs = ['WIDTH', 'FILE_LENGTH', 'X_FIRST', 'X_STEP',
			'Y_FIRST', 'Y_STEP', Z_OFFSET, Z_SCALE, PROJECTION, DATUM]

		for a in attrs:
			self.assertTrue(getattr(self.ras, a) is not None)


	def test_is_dem(self):
		self.assertTrue(hasattr(self.ras, DATUM))
		self.ras = DEM('../../tests/sydney_test/obs/geo_060619-061002.unw')
		self.assertFalse(hasattr(self.ras, DATUM))


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
			_ = self.ras.height_band
			self.fail("Should not be able to access band without open dataset")
		except RasterException:
			pass

		self.ras.open()
		data = self.ras.height_band.ReadAsArray()
		self.assertEqual(data.shape, (72, 47) )


if __name__ == "__main__":
	unittest.main()
