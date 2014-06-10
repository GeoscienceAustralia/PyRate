'''
Library/script to convert ROIPAC headers to ESRI's BIL format.

GDAL lacks a driver to parse ROIPAC headers. This module translates ROIPAC
headers into ESRI's BIL format, which is supported by GDAL. A basic command line
interface is provided for testing purposes. 

The types of ROIPAC files/data used in PyRate are:
* Interferograms: a .unw 32 bit float data file, with a .rsc resource/header.
The binary data is assumed to contain 2 bands, amplitude and phase.

* DEM: with a .unw 16 bit signed int binary data file, and a .rsc header
There is only a single height band for the binary data.

* TODO: describe incidence files, and any others.   


There may be differences with the .rsc file content, with short and long forms.
The short form has 7 fields, covering raster size, location and wavelength. The
longer form can have up to 40 fields (see the test data for examples). PyRate
attempts to handle both forms of header. 

Created on 12/09/2012
@author: Ben Davies, NCI
         ben.davies@anu.edu.au
'''

import os, re
import struct
import datetime

import ifgconstants as ifc
import numpy as np
import gdal, osr
from ifgconstants import PYRATE_DATUM



# ROIPAC RSC header file constants
WIDTH = "WIDTH"
FILE_LENGTH = "FILE_LENGTH"
XMIN = "XMIN"
XMAX = "XMAX"
YMIN = "YMIN"
YMAX = "YMAX"
X_FIRST = "X_FIRST"
X_STEP = "X_STEP"
X_UNIT = "X_UNIT"
Y_FIRST = "Y_FIRST"
Y_STEP = "Y_STEP"
Y_UNIT = "Y_UNIT"
TIME_SPAN_YEAR = "TIME_SPAN_YEAR"

# Old ROIPAC headers (may not be needed)
ORBIT_NUMBER = "ORBIT_NUMBER"
VELOCITY = "VELOCITY"
HEIGHT = "HEIGHT"
EARTH_RADIUS = "EARTH_RADIUS"
WAVELENGTH = "WAVELENGTH"
DATE = "DATE"
DATE12 = "DATE12"
HEADING_DEG = "HEADING_DEG"

# DEM specific
Z_OFFSET = "Z_OFFSET"
Z_SCALE = "Z_SCALE"
PROJECTION = "PROJECTION"
DATUM = "DATUM"

# custom header aliases
MASTER = "MASTER"
SLAVE = "SLAVE"
X_LAST = "X_LAST"
Y_LAST = "Y_LAST"


# store type for each of the header items
INT_HEADERS = [WIDTH, FILE_LENGTH, XMIN, XMAX, YMIN, YMAX, Z_OFFSET, Z_SCALE ]
STR_HEADERS = [X_UNIT, Y_UNIT, ORBIT_NUMBER, DATUM, PROJECTION ]
FLOAT_HEADERS = [X_FIRST, X_STEP, Y_FIRST, Y_STEP, TIME_SPAN_YEAR,
				VELOCITY, HEIGHT, EARTH_RADIUS, WAVELENGTH, HEADING_DEG ]
DATE_HEADERS = [DATE, DATE12]

ROIPAC_HEADER_LEFT_JUSTIFY = 18
ROI_PAC_HEADER_FILE_EXT = "rsc"


# TODO: check for mismatching X,Y cell resolution?
def to_geotiff(hdr, data_path, dest, nodata):
	'Converts GAMMA format data to GeoTIFF image with PyRate metadata'
	is_ifg = ifc.PYRATE_WAVELENGTH_METRES in hdr
	ncols = hdr[ifc.PYRATE_NCOLS]
	nrows = hdr[ifc.PYRATE_NROWS]
	_check_raw_data(is_ifg, data_path, ncols, nrows)
	_check_step_mismatch(hdr)
	
	driver = gdal.GetDriverByName("GTiff")
	dtype = gdal.GDT_Float32 if is_ifg else gdal.GDT_Int16
	ds = driver.Create(dest, ncols, nrows, 1, dtype)
	
	# write custom headers to interferograms
	if is_ifg:
		for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
					ifc.PYRATE_DATE, ifc.PYRATE_DATE2]:
			ds.SetMetadataItem(k, str(hdr[k]))
	
	# position and projection data	
	ds.SetGeoTransform([hdr[ifc.PYRATE_LONG], hdr[ifc.PYRATE_X_STEP], 0,
						hdr[ifc.PYRATE_LAT], 0, hdr[ifc.PYRATE_Y_STEP]])
	
	srs = osr.SpatialReference()
	res = srs.SetWellKnownGeogCS(hdr[ifc.PYRATE_DATUM])
	if res:
		msg = 'Unrecognised projection: %s' % hdr[ifc.PYRATE_DATUM]
		raise RoipacException(msg)
	
	ds.SetProjection(srs.ExportToWkt())
	
	# copy data from the binary file
	band = ds.GetRasterBand(1)
	band.SetNoDataValue(nodata)
		
	if is_ifg:
		fmtstr = '<' + ('f' * ncols) # ifgs are little endian float32s
		bytes_per_col = 4
	else:
		fmtstr = '<' + ('h' * ncols) # DEM is little endian signed int16
		bytes_per_col = 2 
	
	row_bytes = ncols * bytes_per_col

	with open(data_path, 'rb') as f:
		for y in xrange(nrows):
			if is_ifg:
				f.seek(row_bytes, 1) # skip interleaved band 1

			data = struct.unpack(fmtstr, f.read(row_bytes))
			band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)

	ds = None


def _check_raw_data(is_ifg, data_path, ncols, nrows):
	base_size = ncols * nrows
	if is_ifg:
		size = 4 * base_size * 2 # 2 bands of 4 bytes each
	else:
		size = 2 * base_size # single 2 byte band
	
	act_size = os.stat(data_path).st_size 
	if act_size != size:
		msg = '%s should have size %s, not %s. Is the correct file being used?'
		raise RoipacException(msg % (data_path, size, act_size))

def _check_step_mismatch(hdr):
	xs, ys = [abs(i) for i in [hdr[ifc.PYRATE_X_STEP], hdr[ifc.PYRATE_Y_STEP]]] 
	
	if xs != ys:
		msg = 'X and Y cell sizes do not match: %s & %s'
		raise RoipacException(msg % (xs, ys)) 

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
		is_dem = DATUM in headers and 'Z_SCALE' in headers
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
			pass # ignore other headers

	# grab a subset for GeoTIFF conversion
	subset = {}
	subset[ifc.PYRATE_NCOLS] = headers[WIDTH]
	subset[ifc.PYRATE_NROWS] = headers[FILE_LENGTH]
	subset[ifc.PYRATE_LAT] = headers[Y_FIRST]
	subset[ifc.PYRATE_LONG] = headers[X_FIRST]
	subset[ifc.PYRATE_X_STEP] = headers[X_STEP]
	subset[ifc.PYRATE_Y_STEP] = headers[Y_STEP]

	if is_dem:
		subset[ifc.PYRATE_DATUM] = headers[DATUM]
	else:
		subset[ifc.PYRATE_WAVELENGTH_METRES] = headers[WAVELENGTH]

		# grab master/slave dates from header, or the filename		
		has_dates = True if DATE in headers and DATE12 in headers else False
		dates = headers[DATE12] if has_dates else _parse_dates_from(hdr_file)
		subset[ifc.PYRATE_DATE], subset[ifc.PYRATE_DATE2] = dates

		# replace timespan as ROIPAC is ~4 hours different to (slave - master)
		timespan = (subset[ifc.PYRATE_DATE2] - subset[ifc.PYRATE_DATE]).days / ifc.DAYS_PER_YEAR
		subset[ifc.PYRATE_TIME_SPAN] = timespan
		
	# add custom X|Y_LAST for convenience
	subset[X_LAST] = headers[X_FIRST] + (headers[X_STEP] * (headers[WIDTH]))
	subset[Y_LAST] = headers[Y_FIRST] + (headers[Y_STEP] * (headers[FILE_LENGTH]))

	return subset


def _parse_dates_from(filename):
	# process dates from filename if rsc file doesn't have them (skip for DEMs)
	p = re.compile(r'\d{6}-\d{6}') # match 2 sets of 6 digits separated by '-'
	m = p.search(filename)
	
	if m:
		s = m.group()
		min_date_len = 13 # assumes "nnnnnn-nnnnnn" format
		if len(s) == min_date_len:
			return parse_date(s)			
	else:
		msg = "Filename does not include master/slave dates: %s"
		raise RoipacException(msg % filename)


class RoipacException(Exception):
	pass


if __name__ == '__main__':
	import sys
	from optparse import OptionParser
	
	usage = "Usage: %prog [options] ROIPAC_FILE [ROIPAC_FILE...]" 
	parser = OptionParser(usage=usage)

	proj_help = 'GDAL well known projection (eg. "WGS84")'
	parser.add_option('-p', '--projection', help=proj_help, type='str')
	res_help = 'Resource/header file with projection data (usually DEM header)'
	parser.add_option('-r', '--resource-header', help=res_help, metavar='FILE')
	parser.add_option('-n', '--nodata', help='NODATA value', type='float', default=0.0)
	options, args = parser.parse_args()

	if len(args) == 0:
		parser.error(usage)
	
	if options.projection:
		proj = options.projection
	elif options.resource_header:
		hdr = parse_header(options.resource_header)
		if PYRATE_DATUM not in hdr:
			sys.exit('Error: header/resource file does not include DATUM')
		proj = hdr[PYRATE_DATUM]
	else:
		parser.error('Need source for DATUM')

	# translate files
	for path in args:
		try:
			hpath = "%s.%s" % (path, ROI_PAC_HEADER_FILE_EXT)
			hdr = parse_header(hpath)

			if PYRATE_DATUM in hdr:
				raise RoipacException('Datum field already in the header')
			else:
				hdr[PYRATE_DATUM] = proj
			
			dest = "%s.tif" % os.path.splitext(os.path.basename(path))[0]
			to_geotiff(hdr, path, dest, options.nodata)
		except Exception, ex:
			sys.exit(ex.message)			
