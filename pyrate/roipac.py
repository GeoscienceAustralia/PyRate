"""
Utilities for converting ROIPAC headers to ESRI's BIL format.

GDAL lacks a driver to parse ROIPAC headers. This module translates ROIPAC
headers into ESRI's BIL format, which is supported by GDAL. A basic command line
interface is provided for testing purposes.

The types of ROIPAC files/data used in PyRate are:

- *Interferograms*: a .unw 32 bit float data file, with a .rsc resource/header.
  The binary data is assumed to contain 2 bands, amplitude and phase.

- *DEM*: with a .unw 16 bit signed int binary data file, and a .rsc header
  There is only a single height band for the binary data.

.. todo:: Implement & describe incidence files, and any others (for later version)


There may be differences with the .rsc file content, with short and long forms.
The short form has 7 fields, covering raster size, location and wavelength. The
longer form can have up to 40 fields (see the test data for examples). PyRate
attempts to handle both forms of header.

Created on 12/09/2012

.. codeauthor:: Ben Davies, NCI <ben.davies@anu.edu.au>
"""

import os, re, struct, datetime, osr, gdal, sys
import numpy as np
import pyrate.ifgconstants as ifc


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
INT_HEADERS = [WIDTH, FILE_LENGTH, XMIN, XMAX, YMIN, YMAX, Z_OFFSET, Z_SCALE]
STR_HEADERS = [X_UNIT, Y_UNIT, ORBIT_NUMBER, DATUM, PROJECTION ]
FLOAT_HEADERS = [X_FIRST, X_STEP, Y_FIRST, Y_STEP, TIME_SPAN_YEAR,
                VELOCITY, HEIGHT, EARTH_RADIUS, WAVELENGTH, HEADING_DEG]
DATE_HEADERS = [DATE, DATE12]

ROIPAC_HEADER_LEFT_JUSTIFY = 18
ROI_PAC_HEADER_FILE_EXT = "rsc"


def to_geotiff(header, data_path, dest, nodata):
    """Converts GAMMA format data to GeoTIFF image with PyRate metadata"""
    is_ifg = ifc.PYRATE_WAVELENGTH_METRES in header
    ncols = header[ifc.PYRATE_NCOLS]
    nrows = header[ifc.PYRATE_NROWS]
    _check_raw_data(is_ifg, data_path, ncols, nrows)
    _check_step_mismatch(header)

    driver = gdal.GetDriverByName("GTiff")
    dtype = gdal.GDT_Float32 if is_ifg else gdal.GDT_Int16
    ds = driver.Create(dest, ncols, nrows, 1, dtype)

    # write custom headers to interferograms
    if is_ifg:
        for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                    ifc.PYRATE_DATE, ifc.PYRATE_DATE2]:
            ds.SetMetadataItem(k, str(header[k]))

    # position and projection data
    ds.SetGeoTransform([header[ifc.PYRATE_LONG], header[ifc.PYRATE_X_STEP], 0,
                        header[ifc.PYRATE_LAT], 0, header[ifc.PYRATE_Y_STEP]])

    srs = osr.SpatialReference()
    res = srs.SetWellKnownGeogCS(header[ifc.PYRATE_DATUM])
    if res:
        msg = 'Unrecognised projection: %s' % header[ifc.PYRATE_DATUM]
        raise RoipacException(msg)

    ds.SetProjection(srs.ExportToWkt())

    # copy data from the binary file
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    if is_ifg:
        fmtstr = '<' + ('f' * ncols)  # ifgs are little endian float32s
        bytes_per_col = 4
    else:
        fmtstr = '<' + ('h' * ncols)  # DEM is little endian signed int16
        bytes_per_col = 2

    row_bytes = ncols * bytes_per_col

    with open(data_path, 'rb') as f:
        for y in xrange(nrows):
            if is_ifg:
                f.seek(row_bytes, 1)  # skip interleaved band 1

            data = struct.unpack(fmtstr, f.read(row_bytes))
            band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)

    ds = None  # manual close
    del ds


def _check_raw_data(is_ifg, data_path, ncols, nrows):
    base_size = ncols * nrows
    if is_ifg:
        size = 4 * base_size * 2  # 2 bands of 4 bytes each
    else:
        size = 2 * base_size  # single 2 byte band

    act_size = os.stat(data_path).st_size
    if act_size != size:
        msg = '%s should have size %s, not %s. Is the correct file being used?'
        raise RoipacException(msg % (data_path, size, act_size))


def _check_step_mismatch(header):
    xs, ys = [abs(i) for i in [header[ifc.PYRATE_X_STEP], header[ifc.PYRATE_Y_STEP]]]

    if xs != ys:
        msg = 'X and Y cell sizes do not match: %s & %s'
        raise RoipacException(msg % (xs, ys))


def parse_date(dstr):
    """Parses ROI_PAC 'yymmdd' or 'yymmdd-yymmdd' to date or date tuple"""
    def to_date(ds):
        year, month, day = [int(ds[i:i+2]) for i in range(0, 6, 2)]
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
        is_dem = True if DATUM in headers or Z_SCALE in headers or PROJECTION in headers else False
        if is_dem and DATUM not in headers:
            sys.exit('Error: no DATUM parameter in DEM header/resource file')
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
    subset = {ifc.PYRATE_NCOLS: headers[WIDTH],
                ifc.PYRATE_NROWS: headers[FILE_LENGTH],
                ifc.PYRATE_LAT: headers[Y_FIRST],
                ifc.PYRATE_LONG: headers[X_FIRST],
                ifc.PYRATE_X_STEP: headers[X_STEP],
                ifc.PYRATE_Y_STEP: headers[Y_STEP], }

    if is_dem:
        subset[ifc.PYRATE_DATUM] = headers[DATUM]
    else:
        subset[ifc.PYRATE_WAVELENGTH_METRES] = headers[WAVELENGTH]

        # grab master/slave dates from header, or the filename
        has_dates = True if DATE in headers and DATE12 in headers else False
        dates = headers[DATE12] if has_dates else _parse_dates_from(hdr_file)
        subset[ifc.PYRATE_DATE], subset[ifc.PYRATE_DATE2] = dates

        # replace time span as ROIPAC is ~4 hours different to (slave - master)
        timespan = (subset[ifc.PYRATE_DATE2] - subset[ifc.PYRATE_DATE]).days / ifc.DAYS_PER_YEAR
        subset[ifc.PYRATE_TIME_SPAN] = timespan

    # add custom X|Y_LAST for convenience
    subset[X_LAST] = headers[X_FIRST] + (headers[X_STEP] * (headers[WIDTH]))
    subset[Y_LAST] = headers[Y_FIRST] + (headers[Y_STEP] * (headers[FILE_LENGTH]))

    return subset


def _parse_dates_from(filename):
    # process dates from filename if rsc file doesn't have them (skip for DEMs)
    p = re.compile(r'\d{6}-\d{6}')  # match 2 sets of 6 digits separated by '-'
    m = p.search(filename)

    if m:
        s = m.group()
        min_date_len = 13  # assumes "nnnnnn-nnnnnn" format
        if len(s) == min_date_len:
            return parse_date(s)
    else:
        msg = "Filename does not include master/slave dates: %s"
        raise RoipacException(msg % filename)


class RoipacException(Exception):
    pass

