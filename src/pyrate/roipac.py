'''
Tool to convert rasters from ROIPAC to another data format

Created on 12/09/2012
@author: Ben Davies, ANUSF
         ben.davies@anu.edu.au
'''

import os, re, datetime

from ifgconstants import INT_HEADERS, STR_HEADERS, FLOAT_HEADERS, DATE_HEADERS
from ifgconstants import ROI_PAC_HEADER_FILE_EXT, DATE, DATE12, DATUM
from ifgconstants import X_STEP, Y_STEP, FILE_LENGTH, TIME_SPAN_YEAR
from ifgconstants import X_FIRST, X_LAST, Y_FIRST, Y_LAST, WIDTH, MASTER, SLAVE

# constants
ROIPAC_HEADER_LEFT_JUSTIFY = 18
PIXELTYPE_INT = "signedint"
PIXELTYPE_FLOAT = "float"


def filename_pair(base):
	"""Returns tuple of paths: (roi_pac data, roi_pac header file)"""
	b = base.strip()
	return (b, "%s.%s" % (b, ROI_PAC_HEADER_FILE_EXT))


def parse_date(dstr):
	"""Parses ROI_PAC 'yymmdd' or 'yymmdd-yymmdd' to date or date tuple"""
	def to_date(ds):
		year, month, day = [int(ds[i:i+2]) for i in range(0,6,2)]
		year += 1900 if (year <= 99 and year >= 50) else 2000
		return datetime.date(year, month, day)

	if "-" in dstr: # ranged date
		return tuple([to_date(d) for d in dstr.split("-")])
	else:
		return to_date(dstr)


def parse_header(hdr_file):
	"""Parses ROI_PAC header file to a dict"""
	with open(hdr_file) as f:
		text = f.read()

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
		elif k in FLOAT_HEADERS:
			headers[k] = float(headers[k])
		elif k in DATE_HEADERS:
			headers[k] = parse_date(headers[k])
		else:
			raise RoipacException("Unrecognised header element %s: %s " % (k, headers[k]) )
	
	# process dates from filename if rsc file doesn't have them (skip this for DEMs)
	if not headers.has_key(DATUM):
		if headers.has_key(DATE) is False or headers.has_key(DATE12) is False:
			p = re.compile(r'[0-9]+-[0-9]+')
			m = p.search(hdr_file)
			
			if m:
				s = m.group()
				min_date_len = 13 # assumes "nnnnnn-nnnnnn" format
				if len(s) >= min_date_len:
					date12 = parse_date(s)
					headers[DATE] = date12[0]
					headers[DATE12] = date12
			else:
				raise RoipacException("Filename does not include master/slave dates: %s" % hdr_file)

		# add master and slave alias headers
		headers[MASTER] = headers[DATE]
		headers[SLAVE] = headers[DATE12][-1]

		# replace timespan as ROI_PAC is ~4 hours different to (slave - master)
		headers[TIME_SPAN_YEAR] = (headers[SLAVE] - headers[MASTER]).days / 365.25
	
	# add custom X|Y_LAST for convenience
	if not headers.has_key(X_LAST):
		headers[X_LAST] = headers[X_FIRST] + (headers[X_STEP] * (headers[WIDTH]))
	if not headers.has_key(Y_LAST):
		headers[Y_LAST] = headers[Y_FIRST] + (headers[Y_STEP] * (headers[FILE_LENGTH]))
	
	return headers


def to_ehdr_header(hdr, dest=None):
	"""Converts a ROI_PAC header to equivalent EHdr format, allowing GDAL to read
	ROIPAC. 'hdr' is the .rsc header path. 'dest' is alternate path to save to. If
	None, dest is defaulted to the 'base_file_name.hdr' """
	if os.path.isfile(hdr) or os.path.islink(hdr):
		H = parse_header(hdr)
		is_dem = True if H.has_key(DATUM) else False

		if dest is None:
			if is_dem:
				try:
					i = hdr.index("dem.rsc") # assumes ROIPAC has filename.dem & filename.dem.rsc 
					dest = hdr[:i] + "hdr"
				except ValueError as v:
					# DEM possibly not from ROIPAC
					raise NotImplementedError("TODO: handle non ROIPAC filenames?")
			else:
				try:
					i = hdr.index("unw.rsc")
					dest = hdr[:i] + "hdr"
				except ValueError as v:
					raise NotImplementedError("TODO: handle ROIPAC filename errors")

	# calc coords of lower left corner
	yllcorner = H[Y_FIRST] + (H[FILE_LENGTH] * H[Y_STEP])
	if yllcorner > 90 or yllcorner < -90:
		raise RoipacException("Invalid Y latitude for yllcorner: %s" % yllcorner)
	
	# create ESRI/EHdr format header, using ROIPAC defaults
	# NB: ROIPAC uses 0 for phase NODATA, which isn't quite correct. Use zero for
	# now, which allows GDAL to recognise NODATA cells
	with open(dest, "w") as f:
		f.write("ncols %s\n" % H[WIDTH])
		f.write("nrows %s\n" % H[FILE_LENGTH])
		
		# handle cells with different dimensions (square & non-square)
		if H[X_STEP] == abs(H[Y_STEP]):
			f.write("cellsize %s\n" % H[X_STEP])
		else:
			f.write("xdim %s\n" % H[X_STEP])
			f.write("ydim %s\n" % abs(H[Y_STEP]) ) # NB: GDAL reads all zeros if ydim is -ve 		

		f.write("xllcorner %s\n" % H[X_FIRST])
		f.write("yllcorner %s\n" % yllcorner)
		
		if not is_dem:
			f.write("nodata 0\n")
			f.write("layout bil\n") # 1 band DEM doesn't interleave data
		
		f.write("nbands %s\n" % (1 if is_dem else 2) )  # number of bands
		f.write("byteorder lsb\n")
		
		# ROIPAC DEMs are 16 bit signed ints, phase layers are 32 bit floats
		f.write("nbits %s\n" % (16 if is_dem else 32) )
		f.write("pixeltype %s\n" % (PIXELTYPE_INT if is_dem else PIXELTYPE_FLOAT) )

	return dest


def write_roipac_header(params, dest_path):
	"""Writes ROIPAC format header given a dict of parameters"""
	with open(dest_path, 'w') as dest:
		for i in params.items():
			line = i[0].ljust(ROIPAC_HEADER_LEFT_JUSTIFY) + str(i[1]) + "\n"
			dest.write(line)	


class RoipacException(Exception):
	pass
