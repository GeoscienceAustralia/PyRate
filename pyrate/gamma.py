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

import datetime


# constants
GAMMA_DATE = 'date'
GAMMA_FREQUENCY = 'radar_frequency'
SPEED_OF_LIGHT_METRES_PER_SECOND = 3e8


def parse_header(path):
	with open(path) as f:
		text = f.read().splitlines()
		raw_segs = [line.split() for line in text if ':' in line]

	# convert the content into a giant dict of all key, values 
	return dict( (i[0][:-1], i[1:]) for i in raw_segs)

def parse_epoch_header(path):
	'TODO'
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
	'TODO'
	lookup = parse_header(path)
	raise NotImplementedError



def frequency_to_wavelength(freq):
	return SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(path0, path1):
	'Combines'
	
	raise NotImplementedError



class GammaError(Exception):
	pass