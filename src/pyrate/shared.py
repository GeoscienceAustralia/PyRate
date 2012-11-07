'''
Contains objects common to multiple parts of PyRate

Created on 12/09/2012
@author: bpd900
'''

import os
import numpy

try:
	from osgeo import gdal
	from gdalconst import GA_Update
except:
	import gdal, gdalconst

gdal.UseExceptions()

import roipac
from ifgconstants import DATUM

# TODO: add phase_data and amplitude_data properties?
#     Problem: property links to FULL dataset, which may be slower than row by row access
#         row by row access would be efficient, but needes wavelength converted layers



class Ifg(object):
	"""Interferogram class, representing the difference between two acquisitions.
	Ifg objects double as a container for related data."""

	def __init__(self, path, hdr_path=None):
		if hdr_path:
			# handle non default header (eg. for look files with different naming)
			self.data_path, self.hdr_path = path, hdr_path
		else:
			# default the header path
			self.data_path, self.hdr_path = roipac.filename_pair(path)

		# dynamically include header items as class attrs
		header = roipac.parse_header(self.hdr_path)
		self.__dict__.update(header)

		self.ehdr_path = None # path to EHdr format header
		self.dataset = None # for GDAL dataset obj
		self._readonly = None
		self._amp_band = None
		self._phase_band = None

		# TODO: what are these for?
		self.max_variance = None
		self.alpha = None
		self._nan_fraction = None


	def __str__(self):
		return "Ifg('%s')" % self.data_path


	def __repr__(self):
		return "Ifg('%s', '%s')" % (self.data_path, self.hdr_path)


	def open(self, readonly=True):
		'''Opens a interferogram dataset for reading. Creates ESRI/EHdr format
		header in the data dir, so GDAL has access to recognised header.'''
		if self.ehdr_path is None:
			self.ehdr_path = roipac.to_ehdr_header(self.hdr_path)
			args = (self.data_path,) if readonly else (self.data_path, GA_Update)
			self.dataset = gdal.Open(*args)

			if self.dataset is None:
				raise IfgException("Error opening %s" % self.data_path)

		else:
			if self.dataset is not None:
				msg = "open() already called for %s" % self
				raise IfgException(msg)

		self._readonly = readonly
		self.num_cells = self.dataset.RasterYSize * self.dataset.RasterXSize


	@property
	def amp_band(self):
		'''Returns a GDAL Band object for the amplitude band'''
		if self._amp_band is not None:
			return self._amp_band
		else:
			if self.dataset is not None:
				self._amp_band = self.dataset.GetRasterBand(1)
				return self._amp_band
			else:
				raise IfgException("Ifg %s has not been opened" % self.data_path)


	@property
	def phase_band(self):
		'''Returns a GDAL Band object for the phase band'''
		if self._phase_band is not None:
			return self._phase_band
		else:
			if self.dataset is not None:
				self._phase_band = self.dataset.GetRasterBand(2)
				return self._phase_band
			else:
				raise IfgException("Ifg %s has not been opened" % self.data_path)


	@property
	def nan_fraction(self):
		'''Returns 0-1 proportion of NaN cells for the phase band'''

		# TODO: cache nan_count for readonly datasets? Perf benefit vs temp changes to data?
		data = self.phase_band.ReadAsArray()
		nan_count = numpy.sum(numpy.isnan(data))
		return nan_count / float(self.num_cells)



class Raster(object):
	"""Generic raster class for DEMs, initial models etc"""

	def __init__(self, path, hdr_path=None):
		if hdr_path:
			# handle non default header (eg. for look files with different naming)
			self.data_path, self.hdr_path = path, hdr_path
		else:
			# default the header path
			self.data_path, self.hdr_path = roipac.filename_pair(path)

		# dynamically include header items as class attrs
		header = roipac.parse_header(self.hdr_path)
		self.__dict__.update(header)

		self.ehdr_path = None # path to EHdr format header
		self.dataset = None # for GDAL dataset obj
		self._band = None


	def __str__(self):
		return "Raster('%s')" % self.data_path


	def __repr__(self):
		return "Raster('%s', '%s')" % (self.data_path, self.hdr_path)


	def open(self, readonly=True):
		'''Opens generic raster dataset. Creates ESRI/EHdr format header in the data
		dir, so GDAL has a recogniseable header.'''
		if self.ehdr_path is None:
			self.ehdr_path = roipac.to_ehdr_header(self.hdr_path)
			args = (self.data_path,) if readonly else (self.data_path, GA_Update)
			self.dataset = gdal.Open(*args)

			if self.dataset is None:
				raise RasterException("Error opening %s" % self.data_path)

		else:
			if self.dataset is not None:
				msg = "open() already called for %s" % self
				raise RasterException(msg)


	def is_dem(self):
		if hasattr(self, DATUM):
			return True
		return False


	@property
	def band(self):
		if self._band is not None:
			return self._band
		else:
			if self.dataset is not None:
				self._band = self.dataset.GetRasterBand(1)
				return self._band
			else:
				raise RasterException("Ifg %s has not been opened" % self.data_path)



class IfgException(Exception):
	'''Generic exception class for interferogram errors'''
	pass

class RasterException(Exception):
	'''Generic exception for raster errors'''
	pass

class PyRateException(Exception):
	'''Generic exception class for PyRate S/W errors'''
	pass

