'''
Tool to convert rasters from ROIPAC to another data format  

Created on 12/09/2012
@author: Ben Davies, ANUSF
				 ben.davies@anu.edu.au
'''

import os
from shared import IfgConstants


def convert_roipac(src, dest, fmt):
	
	# Two options for conversion:
	# 1) Fake an EHdr header and then read using GDAL normally
	# - is a hack
	# - misses the other header items, these have to be copied in as metadata later, replicating part of #2 anyway  
	#
	# 2) Read header and binary directly
	# - avoids faking the header/cleaner
	#
	# header file has to be read anyway
	
	raise NotImplementedError



def _read_roipac_header(hdr):
	"""Parses ROI_PAC header file into a dict"""
	if os.path.isfile(hdr):
		with open(hdr) as f:
			text = f.read()		
	else:
		text = hdr
	
	lines = [e.split() for e in text.split("\n") if e != ""]
	headers = dict(lines)
	
	for k in headers.keys():
		if k in IfgConstants.INT_HEADERS:
			headers[k] = int(headers[k])
		elif k in IfgConstants.STR_HEADERS:
			headers[k] = str(headers[k])
		else:
			try:
				headers[k] = float(headers[k])
			except ValueError:
				raise RoipacException("Unrecognised floating point header element %s: %s " % (k, headers[k]) )
	
	return headers


class RoipacException(Exception):
	pass 