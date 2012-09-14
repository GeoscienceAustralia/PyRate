'''
Tool to convert rasters from ROIPAC to another data format  

Created on 12/09/2012
@author: Ben Davies, ANUSF
				 ben.davies@anu.edu.au
'''

import os
from osgeo import gdal
from ifgconstants import INT_HEADERS, STR_HEADERS, X_STEP, Y_STEP
from ifgconstants import X_FIRST, Y_FIRST, WIDTH, FILE_LENGTH
from shared import ROI_PAC_HEADER_FILE_EXT



def convert_roipac(src, dest, fmt):
	
	# Two options for conversion:
	# 1) Fake an EHdr header and then read using GDAL normally
	# - is a bit of a hack
	# - misses the other header items, these have to be copied in as metadata later, replicating part of #2 anyway
	# - reuses GDAL for the double band data
	#
	# 2) Read header and binary directly
	# - avoids faking the header/cleaner
	#
	# header file has to be read anyway
	#
	# TODO: read 1 and 2 band files
	
	raise NotImplementedError


def filename_pair(base):
	"""Returns tuple of paths: (roi_pac data, roi_pac header file)"""
	return (base, "%s.%s" % (base, ROI_PAC_HEADER_FILE_EXT))


def _read_roipac_header(hdr):
	"""Parses ROI_PAC header file to a dict"""
	if os.path.isfile(hdr):
		with open(hdr) as f:
			text = f.read()		
	else:
		text = hdr
	
	try:
		lines = [e.split() for e in text.split("\n") if e != ""]
		headers = dict(lines)
	except ValueError:
		raise RoipacException("Unable to parse header content:\n%s" % text)
	
	for k in headers.keys():
		if k in INT_HEADERS:
			headers[k] = int(headers[k])
		elif k in STR_HEADERS:
			headers[k] = str(headers[k])
		else:
			try:
				headers[k] = float(headers[k])
			except ValueError:
				raise RoipacException("Unrecognised header element %s: %s " % (k, headers[k]) )
	
	return headers


def roipac_to_ehdr_header(hdr, dest):
	"""Convenience function to convert a ROI_PAC header to EHdr format. 'hdr' can be
	a path to a header file, or a dict of header elements. 'dest' is path to save to"""
	if os.path.isfile(hdr):
		hdr = _read_roipac_header(hdr)
	
	cellsize = hdr[X_STEP] 
	if cellsize != abs(hdr[Y_STEP]):
		raise RoipacException("Unequal X and Y axis cell sizes: %s, %s" %
												(cellsize, hdr[Y_STEP]) )
	
	# calc coords of lower left corner (EHdr format uses this) 
	yllcorner = hdr[Y_FIRST] + (hdr[FILE_LENGTH] * hdr[Y_STEP])
	if yllcorner > 90 or yllcorner < -90:
		raise RoipacException("Invalid Y latitude for yllcorner: %s" % yllcorner)
		
	with open(dest, "w") as f:
		f.write("ncols %s\n" % hdr[WIDTH])
		f.write("nrows %s\n" % hdr[FILE_LENGTH])
		f.write("cellsize %s\n" % hdr[X_STEP])
		f.write("xllcorner %s\n" % hdr[X_FIRST])
		f.write("yllcorner %s\n" % yllcorner)
		f.write("nbands 2\n")
		f.write("byteorder lsb\n")
		f.write("layout bil\n")
		f.write("nbits 32\n")
		f.write("pixeltype float\n")



class RoipacException(Exception):
	pass
