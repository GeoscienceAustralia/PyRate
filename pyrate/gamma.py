'''
Library/script to convert GAMMA headers to ESRI's BIL format.

GAMMA headers need to be translated into a GDAL recognisable format for use in
PyRate. This module translates GAMMA headers into into ESRI's BIL format,
allowing GDAL to access the raster data

TODO: explain format of GAMMA files (data, DEM, other??)

Created on 29/05/2014
@author: Ben Davies NCI
'''

# TODO: place metadata into GeoTIFF headers with key:values

import struct
import datetime

import gdal, osr
import numpy as np


# constants
GAMMA_DATE = 'date'
GAMMA_FREQUENCY = 'radar_frequency'

GAMMA_WIDTH = 'width'
GAMMA_NROWS = 'nlines'
GAMMA_CORNER_LAT = 'corner_lat'
GAMMA_CORNER_LONG = 'corner_lon'
GAMMA_Y_STEP = 'post_lat'
GAMMA_X_STEP = 'post_lon'

SPEED_OF_LIGHT_METRES_PER_SECOND = 3e8


def to_geotiff(hdr, data_path, dest, nodata):
	'Converts GAMMA format data to GeoTIFF image with PyRate metadata'
	ncols = hdr['NCOLS']
	nrows = hdr['NROWS']
	driver = gdal.GetDriverByName("GTiff") # TODO: hardcoded
	ds = driver.Create(dest, ncols, nrows, 1, gdal.GDT_Float32)
	
	# write custom headers to interferograms
	if hdr.has_key('WAVELENGTH_METRES'):
		for k in hdr:
			ds.SetMetadataItem(k, str(hdr[k]))
	
	# position and projection data	
	gt = [hdr['LONG'], hdr['X_STEP'], 0, hdr['LAT'], 0, hdr['Y_STEP']]
	ds.SetGeoTransform(gt)
	
	# TODO: is this sufficient to geolocate the raster?
	srs = osr.SpatialReference()
	srs.SetWellKnownGeogCS('WGS84')
	ds.SetProjection( srs.ExportToWkt() )	
	
	# copy data from the binary file
	band = ds.GetRasterBand(1)
	band.SetNoDataValue(nodata)
	fmtstr = '!' + ('f' * ncols) # data format is big endian float32s 
	
	with open(data_path, 'rb') as f:
		for y in range(nrows):
			data = struct.unpack(fmtstr, f.read(ncols * 4))
			band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)
	
	ds = None
	

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
	subset['DATE'] = datetime.date(year, month, day)

	# handle conversion to wavelength	
	tmp_freq, unit = lookup[GAMMA_FREQUENCY]
	if unit != "Hz":
		msg = 'Unrecognised unit field for radar_frequency: %s'
		raise GammaError(msg % unit)
		
	subset['WAVELENGTH_METRES'] = frequency_to_wavelength(float(tmp_freq))
	return subset


def parse_dem_header(path):
	'Returns dict of metadata for converting GAMMA to custom PyRate GeoTIFF'
	lookup = parse_header(path)
	subset = {}
	
	# NB: many lookup fields have multiple elements, eg ['1000', 'Hz']
	subset['NCOLS'] = int(lookup[GAMMA_WIDTH][0])
	subset['NROWS'] = int(lookup[GAMMA_NROWS][0])
	
	expected = ['decimal', 'degrees']
	for k in [GAMMA_CORNER_LAT, GAMMA_CORNER_LONG, GAMMA_X_STEP, GAMMA_Y_STEP]:
		units = lookup[GAMMA_CORNER_LAT][1:]
		if  units != expected:
			msg = "Unrecognised units for GAMMA %s field\n. Got %s, expected %s"
			raise GammaError(msg % (k, units, expected))	
	
	subset['LAT'] = float(lookup[GAMMA_CORNER_LAT][0])
	subset['LONG'] = float(lookup[GAMMA_CORNER_LONG][0])
	subset['Y_STEP'] = float(lookup[GAMMA_Y_STEP][0])
	subset['X_STEP'] = float(lookup[GAMMA_X_STEP][0])
	return subset


def frequency_to_wavelength(freq):
	return SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1):
	# TODO: combine dicts from both epoch headers into single ifg header	
	chdr = {}
	
	if hdr1['DATE'] == hdr0['DATE']:
		raise GammaError("Can't combine headers for the same day")
	elif hdr1['DATE'] < hdr0['DATE']:
		raise GammaError("Wrong date order")
		
	chdr['TIME_SPAN_YEAR'] = (hdr1['DATE'] - hdr0['DATE']).days / 365.25
	chdr['DATE'] = hdr0['DATE']
	
	if hdr0['WAVELENGTH_METRES'] != hdr1['WAVELENGTH_METRES']:
		raise GammaError("Wavelengths don't match") 
	
	chdr['WAVELENGTH_METRES'] = hdr0['WAVELENGTH_METRES']  
	return chdr



class GammaError(Exception):
	pass