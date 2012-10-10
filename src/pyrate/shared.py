'''
Created on 12/09/2012

@author: bpd900
'''

import roipac


class Ifg(object):
	"""Interferogram class, representing the difference between two acquisitions
	and other related data."""

	def __init__(self, path):
		self.data_path, self.hdr_path = roipac.filename_pair(path)
		header = roipac.parse_header(self.hdr_path)

		# dynamically include header items as class attributes
		for key, value in header.iteritems():
			if self.__dict__.has_key(key):
				msg = "Attribute %s already exists for %s" % (key, path)
				raise Exception(msg)
			self.__dict__[key] = value

		self.dataset = None # for GDAL dataset

		self.max_variance = None
		self.alpha = None
		self.nodata_fraction = None
