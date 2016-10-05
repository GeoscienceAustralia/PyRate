'''
Tools used for reading GAMMA headers.

.. codeauthor:: Ben Davies, NCI, Sudipta Basak, GA
'''

import os, datetime
import numpy as np
import pyrate.ifgconstants as ifc



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
GAMMA_INCIDENCE = 'incidence_angle'
RADIANS = 'RADIANS'
GAMMA = 'GAMMA'

SPEED_OF_LIGHT_METRES_PER_SECOND = 3e8



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
    """Parses all GAMMA epoch/DEM header file fields into a dictionary"""
    with open(path) as f:
        text = f.read().splitlines()
        raw_segs = [line.split() for line in text if ':' in line]

    # convert the content into a giant dict of all key, values
    return dict((i[0][:-1], i[1:]) for i in raw_segs)


def parse_header_new(path):
    """Parses all GAMMA epoch/DEM header file fields into a dictionary"""
    with open(path) as f:
        text = f.read().splitlines()
        raw_segs = [line.split() for line in text]

    # convert the content into a giant dict of all key, values
    return dict((i[0], i[1]) for i in raw_segs)


def parse_epoch_header(path):
    """Returns dict of the minimum required epoch metadata needed for PyRate"""
    lookup = parse_header(path)
    subset = {}
    year, month, day, hour, mins, sec = [int(float(i))
                                         for i in lookup[GAMMA_DATE][:6]]
    subset[ifc.MASTER_DATE] = datetime.date(year, month, day)
    subset[ifc.MASTER_TIME] = datetime.time(hour, mins, sec)

    # handle conversion of radar frequency to wavelength
    freq, unit = lookup[GAMMA_FREQUENCY]
    if unit != "Hz":
        msg = 'Unrecognised unit field for radar_frequency: %s'
        raise GammaException(msg % unit)

    incidence, unit = lookup[GAMMA_INCIDENCE]
    if unit != "degrees":
        msg = 'Unrecognised unit field for incidence_angle: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_WAVELENGTH_METRES] = frequency_to_wavelength(float(freq))
    subset[ifc.INCIDENCE_ANGLE] = float(incidence)

    return subset


def parse_dem_header(path):
    """Returns dict of metadata for converting GAMMA to custom PyRate GeoTIFF"""
    lookup = parse_header(path)

    # NB: many lookup fields have multiple elements, eg ['1000', 'Hz']
    subset = {ifc.PYRATE_NCOLS: int(lookup[GAMMA_WIDTH][0]),
                ifc.PYRATE_NROWS: int(lookup[GAMMA_NROWS][0])}

    expected = ['decimal', 'degrees']
    for k in [GAMMA_CORNER_LAT, GAMMA_CORNER_LONG, GAMMA_X_STEP, GAMMA_Y_STEP]:
        units = lookup[GAMMA_CORNER_LAT][1:]
        if units != expected:
            msg = "Unrecognised units for GAMMA %s field\n. Got %s, expected %s"
            raise GammaException(msg % (k, units, expected))

    subset[ifc.PYRATE_LAT] = float(lookup[GAMMA_CORNER_LAT][0])
    subset[ifc.PYRATE_LONG] = float(lookup[GAMMA_CORNER_LONG][0])
    subset[ifc.PYRATE_Y_STEP] = float(lookup[GAMMA_Y_STEP][0])
    subset[ifc.PYRATE_X_STEP] = float(lookup[GAMMA_X_STEP][0])
    subset[ifc.PYRATE_DATUM] = "".join(lookup[GAMMA_DATUM])
    subset[ifc.PYRATE_INSAR_PROCESSOR] = GAMMA
    return subset


def frequency_to_wavelength(freq):
    """
    Convert radar frequency to wavelength
    """
    return SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1, dem_hdr):
    """
    Combines both epoch header lookups into single ifg header/dict

    hdr0: header for the earliest/master ifg
    hdr1: header for the latest/slave ifg
    dem_hdr: dict of DEM header attributes
    """
    if not all([isinstance(a, dict) for a in [hdr0, hdr1, dem_hdr]]):
        raise GammaException('Header args need to be dicts')

    date0, date1 = hdr0[ifc.MASTER_DATE], hdr1[ifc.MASTER_DATE]
    if date0 == date1:
        raise GammaException("Can't combine headers for the same day")
    elif date1 < date0:
        raise GammaException("Wrong date order")

    chdr = {ifc.PYRATE_TIME_SPAN: (date1 - date0).days / ifc.DAYS_PER_YEAR,
            ifc.MASTER_DATE: date0,
            ifc.MASTER_TIME: hdr0[ifc.MASTER_TIME],
            ifc.SLAVE_DATE: date1,
            ifc.SLAVE_TIME: hdr1[ifc.MASTER_TIME],
            ifc.PYRATE_PHASE_UNITS: RADIANS,
            ifc.PYRATE_INSAR_PROCESSOR: GAMMA}

    wavelen = hdr0[ifc.PYRATE_WAVELENGTH_METRES]
    if np.isclose(wavelen, hdr1[ifc.PYRATE_WAVELENGTH_METRES], atol=1e-6):
        chdr[ifc.PYRATE_WAVELENGTH_METRES] = wavelen
    else:
        args = (chdr[ifc.MASTER_DATE], chdr[ifc.SLAVE_DATE])
        msg = "Wavelength mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)
    # non-cropped, non-multilooked geotif process step information added
    chdr[ifc.PROCESS_STEP] = ifc.GEOTIFF

    chdr.update(dem_hdr)  # add geographic data
    return chdr


def manage_headers(demHeaderFile, headerPaths):
    demHeader = parse_dem_header(demHeaderFile)
    # find param files containing filename dates
    if len(headerPaths) == 2:
        headers = [parse_epoch_header(hp) for hp in headerPaths]
        combinedHeader = combine_headers(headers[0], headers[1], demHeader)
    else:
        # probably have DEM or incidence file
        combinedHeader = demHeader
        combinedHeader[ifc.PROCESS_STEP] = ifc.GEOTIFF

    return combinedHeader


class GammaException(Exception):
    pass

