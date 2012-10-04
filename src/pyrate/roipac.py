'''
Tool to convert rasters from ROIPAC to another data format

Created on 12/09/2012
@author: Ben Davies, ANUSF
         ben.davies@anu.edu.au
'''

import os, re, datetime

from ifgconstants import INT_HEADERS, STR_HEADERS, FLOAT_HEADERS, DATE_HEADERS
from ifgconstants import ROI_PAC_HEADER_FILE_EXT, DATE, DATE12
from ifgconstants import X_STEP, Y_STEP, FILE_LENGTH, TIME_SPAN_YEAR
from ifgconstants import X_FIRST, Y_FIRST, WIDTH, MASTER, SLAVE



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

	if headers.has_key(DATE) is False or headers.has_key(DATE12) is False:
		# probably have short form header without dates, get date from path
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
	return headers


def to_ehdr_header(hdr, dest=None):
	"""Converts a ROI_PAC header to equivalent EHdr format, allowing GDAL to read
	ROIPAC. 'hdr' is the .rsc header path. 'dest' is an alternate path to save to"""
	if os.path.isfile(hdr) or os.path.islink(hdr):
		H = parse_header(hdr)
		if dest is None:
			i = hdr.index("unw.rsc")
			dest = hdr[:i] + "hdr"

	cellsize = H[X_STEP]
	if cellsize != abs(H[Y_STEP]):
		err = "Unequal X and Y axis cell sizes: %s, %s" % (cellsize, H[Y_STEP])
		raise RoipacException(err)

	# calc coords of lower left corner (EHdr format uses this)
	yllcorner = H[Y_FIRST] + (H[FILE_LENGTH] * H[Y_STEP])
	if yllcorner > 90 or yllcorner < -90:
		raise RoipacException("Invalid Y latitude for yllcorner: %s" % yllcorner)

	with open(dest, "w") as f:
		f.write("ncols %s\n" % H[WIDTH])
		f.write("nrows %s\n" % H[FILE_LENGTH])
		f.write("cellsize %s\n" % H[X_STEP])
		f.write("xllcorner %s\n" % H[X_FIRST])
		f.write("yllcorner %s\n" % yllcorner)
		f.write("nbands 2\n")
		f.write("byteorder lsb\n")
		f.write("layout bil\n")
		f.write("nbits 32\n")
		f.write("pixeltype float\n")


class RoipacException(Exception):
	pass
