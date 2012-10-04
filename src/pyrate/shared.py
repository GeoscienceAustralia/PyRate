'''
Created on 12/09/2012

@author: bpd900
'''

import roipac


class Ifg(object):
	"""TODO:"""

	def __init__(self, path):
		self.data_path, self.hdr_path = roipac.filename_pair(path)
		header = roipac.parse_header(self.hdr_path)

		# dynamically include header items as class attributes
		for key, value in header.iteritems():
			if self.__dict__.has_key(key):
				raise Exception("Attribute %s exists for Interferogram %s" % (key, path))
			self.__dict__[key] = value

		self.dataset = None

		self.max_variance = None
		self.alpha = None
		self.nodata_fraction = None
