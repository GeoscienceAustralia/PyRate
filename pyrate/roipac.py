'''
Library/script to convert ROIPAC headers to ESRI's BIL format.

GDAL cannot parse ROIPAC headers, preventing interoperability. This module
translates ROIPAC headers into ESRI's BIL format, which is supported by GDAL. A
command line interface are provided for testing purposes. 

Created on 12/09/2012
@author: Ben Davies, NCI
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
		msg = "Unable to parse content of %s. Is it a ROIPAC header file?"
		raise RoipacException(msg % hdr_file)

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
			msg = "Unrecognised header element %s: %s "
			raise RoipacException(msg % (k, headers[k]) )

	# process dates from filename if rsc file doesn't have them (skip for DEMs)
	if not headers.has_key(DATUM):
		if headers.has_key(DATE) is False or headers.has_key(DATE12) is False:
			p = re.compile(r'\d{6}-\d{6}') # match 2 sets of 6 digits separated by '-'
			m = p.search(hdr_file)

			if m:
				s = m.group()
				min_date_len = 13 # assumes "nnnnnn-nnnnnn" format
				if len(s) == min_date_len:
					date12 = parse_date(s)
					headers[DATE] = date12[0]
					headers[DATE12] = date12
			else:
				msg = "Filename does not include master/slave dates: %s"
				raise RoipacException(msg % hdr_file)

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
	"""
	Converts ROI_PAC header files to equivalent ESRI BIL/GDAL EHdr format. This
	allows GDAL to recognise and read ROIPAC datasets.
	
	hdr: path to the .rsc header file.
	dest: path to save new header to, if None, defaults to 'base_file_name.hdr'
	"""
	if os.path.isfile(hdr) or os.path.islink(hdr):
		H = parse_header(hdr)
		is_dem = True if H.has_key(DATUM) else False

		if dest is None:
			# determine default destination file
			if is_dem:
				try:
					# assumes ROIPAC uses filename.dem & filename.dem.rsc
					i = hdr.index("dem.rsc")
					dest = hdr[:i] + "hdr"
				except ValueError:
					# DEM probably not from ROIPAC
					msg = "Unrecognised file naming pattern for %s"
					raise RoipacException(msg % hdr)
			else:
				i = max(hdr.rfind("unw.rsc"), hdr.rfind("tif.rsc"))
				if i > 0:
					dest = hdr[:i] + "hdr"
				else:
					msg = "Unrecognised file naming pattern for %s"
					raise RoipacException(msg % hdr)
	else:
		raise IOError("%s not a valid header file" % hdr)

	# calc coords of lower left corner
	yllcorner = H[Y_FIRST] + (H[FILE_LENGTH] * H[Y_STEP])
	if yllcorner > 90 or yllcorner < -90:
		raise RoipacException("Invalid Y latitude for yllcorner: %s" % yllcorner)

	# create ESRI BIL format header, using ROIPAC defaults
	# TODO: ROIPAC uses 0 for phase NODATA, which isn't quite correct. Use zero
	# for now to allow GDAL to recognise NODATA cells
	with open(dest, "w") as f:
		f.write("ncols %s\n" % H[WIDTH])
		f.write("nrows %s\n" % H[FILE_LENGTH])

		# handle cells with different dimensions (square & non-square)
		if H[X_STEP] == abs(H[Y_STEP]):
			f.write("cellsize %s\n" % H[X_STEP])
		else:
			f.write("xdim %s\n" % H[X_STEP])
			# GDAL reads zeros if ydim is negative
			f.write("ydim %s\n" % abs(H[Y_STEP]) )

		f.write("xllcorner %s\n" % H[X_FIRST])
		f.write("yllcorner %s\n" % yllcorner)
		f.write("byteorder lsb\n")

		if not is_dem:
			f.write("nodata 0\n")
			f.write("layout bil\n") # 1 band DEM doesn't interleave data
			f.write("nbands 2\n")  # number of bands

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


if __name__ == '__main__':
	import sys
	usage = "Usage: rp2ehdr.py [ROIPAC unw.rsc files ...]\n" 
	if len(sys.argv) < 2:
		sys.stderr.write(usage)
		sys.exit()

	for path in sys.argv[1:]:
		try:
			to_ehdr_header(path)
		except Exception as ex:
			sys.exit(ex.message)
