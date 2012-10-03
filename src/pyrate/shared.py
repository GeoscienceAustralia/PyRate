'''
Created on 12/09/2012

@author: bpd900
'''

import roipac


class Ifg(object):
	"""TODO:"""

	def __init__(self, path):
		self.data_path, self.hdr_path = roipac.filename_pair(path)
		self.header = roipac.parse_header(self.hdr_path)
		self.dataset = None

		self.max_variance = None
		self.alpha = None
		self.nodata_fraction = None
