#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains tools for reading ROI_PAC format input data.
"""
import os
from pathlib import Path
import re
import datetime
import pyrate.core.ifgconstants as ifc
from pyrate.core import config as cf
from pyrate.core.shared import extract_epochs_from_filename

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
FIRST = "FIRST"
SECOND = "SECOND"
X_LAST = "X_LAST"
Y_LAST = "Y_LAST"
RADIANS = "RADIANS"
ROIPAC = "ROIPAC"

# store type for each of the header items
INT_HEADERS = [WIDTH, FILE_LENGTH, XMIN, XMAX, YMIN, YMAX, Z_OFFSET, Z_SCALE]
STR_HEADERS = [X_UNIT, Y_UNIT, ORBIT_NUMBER, DATUM, PROJECTION]
FLOAT_HEADERS = [X_FIRST, X_STEP, Y_FIRST, Y_STEP, TIME_SPAN_YEAR,
                 VELOCITY, HEIGHT, EARTH_RADIUS, WAVELENGTH, HEADING_DEG]
DATE_HEADERS = [DATE, DATE12]

ROIPAC_HEADER_LEFT_JUSTIFY = 18
ROI_PAC_HEADER_FILE_EXT = ".rsc"

def parse_date(dstr):
    """
    Parses ROI_PAC 'yymmdd' or 'yymmdd-yymmdd' format string to datetime.

    :param str dstr: 'date' or 'date1-date2' string

    :return: dstr: datetime string or tuple
    :rtype: str or tuple
    """
    def to_date(date_str):
        """convert string to datetime"""
        year, month, day = [int(date_str[i:i+2]) for i in range(0, 6, 2)]
        year += 1900 if ((year <= 99) and (year >= 50)) else 2000
        return datetime.date(year, month, day)

    if "-" in dstr:  # ranged date
        return tuple([to_date(d) for d in dstr.split("-")])
    else:
        return to_date(dstr)


def parse_header(hdr_file):
    """
    Parses ROI_PAC header file metadata to a dictionary.

    :param str hdr_file: `path to ROI_PAC *.rsc file`

    :return: subset: subset of metadata
    :rtype: dict
    """
    with open(hdr_file, encoding="utf8", errors='ignore') as f:
        text = f.read()

    try:
        lines = [e.split() for e in text.split("\n") if e != ""]
        headers = dict(lines)
        is_dem = True if DATUM in headers or Z_SCALE in headers \
                         or PROJECTION in headers else False
        if is_dem and DATUM not in headers:
            msg = 'No "DATUM" parameter in DEM header/resource file'
            raise RoipacException(msg)
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
        else:  # pragma: no cover
            pass  # ignore other headers

    # grab a subset for GeoTIFF conversion
    subset = {ifc.PYRATE_NCOLS: headers[WIDTH],
              ifc.PYRATE_NROWS: headers[FILE_LENGTH],
              ifc.PYRATE_LAT: headers[Y_FIRST],
              ifc.PYRATE_LONG: headers[X_FIRST],
              ifc.PYRATE_X_STEP: headers[X_STEP],
              ifc.PYRATE_Y_STEP: headers[Y_STEP]}

    if is_dem:
        subset[ifc.PYRATE_DATUM] = headers[DATUM]
    else:
        subset[ifc.PYRATE_WAVELENGTH_METRES] = headers[WAVELENGTH]

        # grab first/second dates from header, or the filename
        has_dates = True if DATE in headers and DATE12 in headers else False
        dates = headers[DATE12] if has_dates else _parse_dates_from(hdr_file)
        subset[ifc.FIRST_DATE], subset[ifc.SECOND_DATE] = dates

        # replace time span as ROIPAC is ~4 hours different to (second minus first)
        timespan = (subset[ifc.SECOND_DATE] - subset[ifc.FIRST_DATE]).days / ifc.DAYS_PER_YEAR
        subset[ifc.PYRATE_TIME_SPAN] = timespan

        # Add data units of interferogram
        subset[ifc.DATA_UNITS] = RADIANS

    # Add InSAR processor flag
    subset[ifc.PYRATE_INSAR_PROCESSOR] = ROIPAC

    # add custom X|Y_LAST for convenience
    subset[X_LAST] = headers[X_FIRST] + (headers[X_STEP] * (headers[WIDTH]))
    subset[Y_LAST] = headers[Y_FIRST] + (headers[Y_STEP] * (headers[FILE_LENGTH]))

    return subset


def _parse_dates_from(filename):
    """Determine dates from file name"""
    # pylint: disable=invalid-name
    # process dates from filename if rsc file doesn't have them (skip for DEMs)
    p = re.compile(r'\d{6}-\d{6}')  # match 2 sets of 6 digits separated by '-'
    m = p.search(filename)

    if m:
        s = m.group()
        min_date_len = 13  # assumes "nnnnnn-nnnnnn" format
        if len(s) == min_date_len:
            return parse_date(s)
    else:  # pragma: no cover
        msg = "Filename does not include first/second image dates: %s"
        raise RoipacException(msg % filename)


def manage_header(header_file, projection):
    """
    Manage header files for ROI_PAC interferograms and DEM files.
    NB: projection = roipac.parse_header(dem_file)[ifc.PYRATE_DATUM]

    :param str header_file: `ROI_PAC *.rsc header file path`
    :param projection: Projection obtained from dem header.

    :return: combined_header: Combined metadata dictionary
    :rtype: dict
    """
    header = parse_header(header_file)
    if ifc.PYRATE_DATUM not in header:  # DEM already has DATUM
        header[ifc.PYRATE_DATUM] = projection
    header[ifc.DATA_TYPE] = ifc.ORIG  # non-cropped, non-multilooked geotiff
    return header


def roipac_header(file_path, params):
    """
    Function to obtain a header for roipac interferogram file or converted
    geotiff.
    """
    rsc_file = params[cf.DEM_HEADER_FILE]
    p = Path(file_path)
    if rsc_file is not None:
        projection = parse_header(rsc_file)[ifc.PYRATE_DATUM]
    else:
        raise RoipacException('No DEM resource/header file is provided')
    if file_path.endswith('_dem.tif'):
        header_file = os.path.join(params[cf.DEM_HEADER_FILE])
    elif file_path.endswith('unw_ifg.tif') or file_path.endswith('unw.tif'):
        # TODO: improve this
        interferogram_epoches = extract_epochs_from_filename(p.name)
        for header_path in params[cf.HEADER_FILE_PATHS]:
            h = Path(header_path.unwrapped_path)
            header_epochs = extract_epochs_from_filename(h.name)
            if set(header_epochs).__eq__(set(interferogram_epoches)):
                header_file = header_path.unwrapped_path
                break
    else:
        header_file = "%s%s" % (file_path, ROI_PAC_HEADER_FILE_EXT)

    header = manage_header(header_file, projection)

    return header


class RoipacException(Exception):
    """
    Convenience class for throwing exception
    """
