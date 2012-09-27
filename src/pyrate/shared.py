'''
Created on 12/09/2012

@author: bpd900
'''

import roipac



class Ifg(object):
	"""TODO:"""
	
	def __init__(self, path):
		self.data_path, self.header_path = roipac.filename_pair(path)
		self.header = roipac.parse_header(self.header_path)
		self.dataset = None
