'''
Library/script to convert GAMMA headers to ESRI's BIL format.

GAMMA headers need to be translated into a GDAL recognisable format for use in
PyRate. This module translates GAMMA headers into into ESRI's BIL format,
allowing GDAL to access the raster data

The types of GAMMA files converted by PyRate are:
* DEM: with a .unw float32 binary data file (MSB order), & '.par' header. There
is only a single height band in the binary data.

* Interferograms: have a .unw float32 bit binary data file (MSB order), with a
single band for phase data. Two .par resource/header files, each containing
details of the epochs used to create the interferogram. No geographic date is 
sotred in these, so the DEM header is required for raster sizes/location etc.

* TODO: describe incidence files (and any others (for later versions)

Created on 29/05/2014
@author: Ben Davies, NCI
'''

import os
import re
import struct
import datetime
import ifgconstants as ifc 
from glob import glob
from os.path import join

import gdal, osr
import numpy as np


# constants
GAMMA_DATE = 'date'
GAMMA_WIDTH = 'width'
GAMMA_NROWS = 'nlines'
GAMMA_CORNER_LAT = 'corner_lat'
GAMMA_CORNER_LONG = 'corner_lon'
GAMMA_Y_STEP = 'post_lat'
GAMMA_X_STEP = 'post_lon'
GAMMA_DATUM = 'ellipsoid_name'
GAMMA_FREQUENCY = 'radar_frequency'


SPEED_OF_LIGHT_METRES_PER_SECOND = 3e8



def to_geotiff(hdr, data_path, dest, nodata):
	'Converts GAMMA format data to GeoTIFF image with PyRate metadata'
	is_ifg = hdr.has_key(ifc.PYRATE_WAVELENGTH_METRES)
	ncols = hdr[ifc.PYRATE_NCOLS]
	nrows = hdr[ifc.PYRATE_NROWS]
	_check_raw_data(data_path, ncols, nrows)
	_check_step_mismatch(hdr)

	driver = gdal.GetDriverByName("GTiff")
	ds = driver.Create(dest, ncols, nrows, 1, gdal.GDT_Float32)

	# write custom headers to interferograms
	if is_ifg:
		for k in (ifc.PYRATE_DATE, ifc.PYRATE_DATE2,
				ifc.PYRATE_TIME_SPAN, ifc.PYRATE_WAVELENGTH_METRES):
			ds.SetMetadataItem(k, str(hdr[k]))

	# position and projection data
	ds.SetGeoTransform([hdr[ifc.PYRATE_LONG], hdr[ifc.PYRATE_X_STEP], 0,
						hdr[ifc.PYRATE_LAT], 0, hdr[ifc.PYRATE_Y_STEP]])

	srs = osr.SpatialReference()
	res = srs.SetWellKnownGeogCS(hdr[ifc.PYRATE_DATUM])

	if res:
		msg = 'Unrecognised projection: %s' % hdr[ifc.PYRATE_DATUM]
		raise GammaException(msg)

	ds.SetProjection(srs.ExportToWkt())

	# copy data from the binary file
	band = ds.GetRasterBand(1)
	band.SetNoDataValue(nodata)
	fmtstr = '!' + ('f' * ncols) # data format is big endian float32s

	with open(data_path, 'rb') as f:
		for y in range(nrows):
			data = struct.unpack(fmtstr, f.read(ncols * 4))
			band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)

	ds = None

def _check_raw_data(data_path, ncols, nrows):
	size = ncols * nrows * 4 # DEM and Ifg data are 4 byte floats 
	act_size = os.stat(data_path).st_size
	if act_size != size:
		msg = '%s should have size %s, not %s. Is the correct file being used?'
		raise GammaException(msg % (data_path, size, act_size))

def _check_step_mismatch(hdr):
	xs, ys = [abs(i) for i in [hdr[ifc.PYRATE_X_STEP], hdr[ifc.PYRATE_Y_STEP]]] 

	if xs != ys:
		msg = 'X and Y cell sizes do not match: %s & %s'
		raise GammaException(msg % (xs, ys)) 

def parse_header(path):
	'Parses all GAMMA epoch/DEM header file fields into a dictionary'
	with open(path) as f:
		text = f.read().splitlines()
		raw_segs = [line.split() for line in text if ':' in line]

	# convert the content into a giant dict of all key, values
	return dict( (i[0][:-1], i[1:]) for i in raw_segs)

def parse_epoch_header(path):
	'Returns dict of the minimum required epoch metadata needed for PyRate'
	lookup = parse_header(path)

	subset = {}
	year, month, day = [int(i) for i in lookup[GAMMA_DATE]]
	subset[ifc.PYRATE_DATE] = datetime.date(year, month, day)

	# handle conversion to wavelength	
	freq, unit = lookup[GAMMA_FREQUENCY]
	if unit != "Hz":
		msg = 'Unrecognised unit field for radar_frequency: %s'
		raise GammaException(msg % unit)

	subset[ifc.PYRATE_WAVELENGTH_METRES] = frequency_to_wavelength(float(freq))
	return subset


def parse_dem_header(path):
	'Returns dict of metadata for converting GAMMA to custom PyRate GeoTIFF'
	lookup = parse_header(path)
	subset = {}

	# NB: many lookup fields have multiple elements, eg ['1000', 'Hz']
	subset[ifc.PYRATE_NCOLS] = int(lookup[GAMMA_WIDTH][0])
	subset[ifc.PYRATE_NROWS] = int(lookup[GAMMA_NROWS][0])

	expected = ['decimal', 'degrees']
	for k in [GAMMA_CORNER_LAT, GAMMA_CORNER_LONG, GAMMA_X_STEP, GAMMA_Y_STEP]:
		units = lookup[GAMMA_CORNER_LAT][1:]
		if  units != expected:
			msg = "Unrecognised units for GAMMA %s field\n. Got %s, expected %s"
			raise GammaException(msg % (k, units, expected))	

	subset[ifc.PYRATE_LAT] = float(lookup[GAMMA_CORNER_LAT][0])
	subset[ifc.PYRATE_LONG] = float(lookup[GAMMA_CORNER_LONG][0])
	subset[ifc.PYRATE_Y_STEP] = float(lookup[GAMMA_Y_STEP][0])
	subset[ifc.PYRATE_X_STEP] = float(lookup[GAMMA_X_STEP][0])
	subset[ifc.PYRATE_DATUM] = "".join(lookup[GAMMA_DATUM])
	return subset


def frequency_to_wavelength(freq):
	return SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1, dem_hdr):
	'''
	Combines both epoch header lookups into single ifg header/dict
	
	hdr0: header for the earliest/master ifg
	hdr1: header for the latest/slave ifg
	dem_hdr: dict of DEM header attributes
	'''
	if not all([isinstance(a, dict) for a in [hdr0, hdr1, dem_hdr]]):
		raise GammaException('Header args need to be dicts')

	chdr = {}
	date0, date1 = hdr0[ifc.PYRATE_DATE], hdr1[ifc.PYRATE_DATE]
	if date0 == date1:
		raise GammaException("Can't combine headers for the same day")
	elif date1 < date0:
		raise GammaException("Wrong date order")

	chdr[ifc.PYRATE_TIME_SPAN] = (date1 - date0).days / ifc.DAYS_PER_YEAR
	chdr[ifc.PYRATE_DATE] = date0
	chdr[ifc.PYRATE_DATE2] = date1 # add 2nd date as it may not be in file name

	wavelen = hdr0[ifc.PYRATE_WAVELENGTH_METRES] 
	if wavelen == hdr1[ifc.PYRATE_WAVELENGTH_METRES]:
		chdr[ifc.PYRATE_WAVELENGTH_METRES] = wavelen
	else:
		args = (chdr[ifc.PYRATE_DATE], chdr[ifc.PYRATE_DATE2])
		msg = "Wavelength mismatch, check both header files for %s & %s"
		raise GammaException(msg % args)

	chdr.update(dem_hdr) # add geographic data
	return chdr


class GammaException(Exception):
	pass



def main():
	import sys
	from optparse import OptionParser

	usage = 'Usage: %prog [options] DEM-HEADER GAMMA_FILE [GAMMA_FILE...]' 
	parser = OptionParser(usage=usage)
	parser.add_option('-n', '--nodata', help='NODATA value', type='float', default=0.0)
	parser.add_option('-d', '--dest-dir', help='Write output to DIR', type='string')
	options, args = parser.parse_args()

	if len(args) < 2:
		parser.error(usage)

	dem_hdr = parse_dem_header(args[0])

	for path in args[1:]:
		dest = '%s.tif' % os.path.splitext(os.path.basename(path))[0]

		if options.dest_dir:
			dest = join(options.dest_dir, dest)

		if os.path.exists(dest):
			sys.exit('Error: %s already exists' % dest)

		try:
			# find param files contaning filename dates
			ptn = re.compile(r'\d{8}') # match 8 digits for the date
			matches = ptn.findall(path)

			if len(matches) == 2:
				srch_dir = os.path.split(path)[0]
				hpaths = [glob(join(srch_dir, '*%s*.par' % m))[0] for m in matches]
				tmp = [parse_epoch_header(hp) for hp in hpaths]
				hdrs = combine_headers(tmp[0], tmp[1], dem_hdr)
			else:
				# probably have DEM or incidence file
				hdrs = dem_hdr

			to_geotiff(hdrs, path, dest, options.nodata)
		except Exception, ex:
			sys.exit('Error: %s' % ex.message)


if __name__ == '__main__':
	main()
