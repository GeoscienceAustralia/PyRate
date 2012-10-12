'''
Created on 12/09/2012

@author: bpd900
'''

import os
import gdal

import roipac


class Ifg(object):
	"""Interferogram class, representing the difference between two acquisitions
	and other related data."""

	def __init__(self, path):
		self.data_path, self.hdr_path = roipac.filename_pair(path)
		header = roipac.parse_header(self.hdr_path)
		self.ehdr_path = None # path to EHdr format header

		# dynamically include header items as class attributes
		for key, value in header.iteritems():
			if self.__dict__.has_key(key):
				msg = "Attribute %s already exists for %s" % (key, path)
				raise Exception(msg)
			self.__dict__[key] = value

		self.dataset = None # for GDAL dataset obj
		self.band = None

		# TODO: what are these for?
		self.max_variance = None
		self.alpha = None
		self.nodata_fraction = None


	def __str__(self):
		return "Ifg('%s')" % self.data_path


	def open(self):
		'''Opens a interferogram dataset for reading. Creates ESRI/EHdr format
		header in the data dir, so GDAL has access to recognised header.'''
		if self.ehdr_path is None:
			self.ehdr_path = roipac.to_ehdr_header(self.hdr_path)
			self.dataset = gdal.Open(self.data_path)
		else:
			if self.dataset is not None:
				msg = "open() already called for %s" % self
				raise IfgException(msg)


class IfgException(Exception):
	'''Generic exception class for interferogram errors'''
	pass

class PyRateException(Exception):
	'''Generic exception class for PyRate S/W errors'''
	pass

