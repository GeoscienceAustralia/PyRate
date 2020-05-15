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
This Python module contains tools for reading GAMMA format input data.
"""
# coding: utf-8
import re
import os
from pathlib import Path
from datetime import date, time, timedelta
import numpy as np
import pyrate.core.ifgconstants as ifc
from pyrate.core import config as cf
from pyrate.core.shared import extract_epochs_from_filename

# constants
GAMMA_DATE = 'date'
GAMMA_TIME = 'center_time'
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


def _parse_header(path):
    """Parses all GAMMA header file fields into a dictionary"""
    with open(path) as f:
        text = f.read().splitlines()
        raw_segs = [line.split() for line in text if ':' in line]

    # convert the content into a giant dict of all key, values
    return dict((i[0][:-1], i[1:]) for i in raw_segs)


def parse_epoch_header(path):
    """
    Returns dictionary of epoch metadata required for PyRate

    :param str path: `Full path to Gamma *slc.par file`

    :return: subset: subset of full metadata
    :rtype: dict
    """
    lookup = _parse_header(path)
    subset = _parse_date_time(lookup)

    # handle conversion of radar frequency to wavelength
    freq, unit = lookup[GAMMA_FREQUENCY]
    if unit != "Hz":  # pragma: no cover
        msg = 'Unrecognised unit field for radar_frequency: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_WAVELENGTH_METRES] = _frequency_to_wavelength(float(freq))

    incidence, unit = lookup[GAMMA_INCIDENCE]
    if unit != "degrees":  # pragma: no cover
        msg = 'Unrecognised unit field for incidence_angle: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_INCIDENCE_DEGREES] = float(incidence)

    return subset


def _parse_date_time(lookup):
    """Grab date and time metadata and convert to datetime objects"""
    subset = {}
    if len(lookup[GAMMA_DATE]) == 3:  # pragma: no cover
        year, month, day, = [int(float(i)) for i in lookup[GAMMA_DATE][:3]]
        if lookup.get(GAMMA_TIME) is not None:
            t = lookup[GAMMA_TIME][0]
            h, m, s = str(timedelta(seconds=float(t))).split(":")
            hour = int(h)
            min = int(m)
            sec = int(s.split(".")[0])
        else:
            # Occasionally GAMMA header has no time information - default to midnight
            hour, min, sec = 0, 0, 0
    elif len(lookup[GAMMA_DATE]) == 6:
        year, month, day, hour, min, sec = [int(float(i)) for i in lookup[GAMMA_DATE][:6]]
    else:  # pragma: no cover
        msg = "Date and time information not complete in GAMMA headers"
        raise GammaException(msg)

    subset[ifc.MASTER_DATE] = date(year, month, day)
    subset[ifc.MASTER_TIME] = time(hour, min, sec)

    return subset


def parse_dem_header(path):
    """
    Returns dictionary of DEM metadata required for PyRate

    :param str path: `Full path to Gamma *dem.par file`

    :return: subset: subset of full metadata
    :rtype: dict
    """
    lookup = _parse_header(path)

    # NB: many lookup fields have multiple elements, eg ['1000', 'Hz']
    subset = {ifc.PYRATE_NCOLS: int(lookup[GAMMA_WIDTH][0]), ifc.PYRATE_NROWS: int(lookup[GAMMA_NROWS][0])}

    expected = ['decimal', 'degrees']
    for k in [GAMMA_CORNER_LAT, GAMMA_CORNER_LONG, GAMMA_X_STEP, GAMMA_Y_STEP]:
        units = lookup[GAMMA_CORNER_LAT][1:]
        if units != expected:  # pragma: no cover
            msg = "Unrecognised units for GAMMA %s field\n. Got %s, expected %s"
            raise GammaException(msg % (k, units, expected))

    subset[ifc.PYRATE_LAT] = float(lookup[GAMMA_CORNER_LAT][0])
    subset[ifc.PYRATE_LONG] = float(lookup[GAMMA_CORNER_LONG][0])
    subset[ifc.PYRATE_Y_STEP] = float(lookup[GAMMA_Y_STEP][0])
    subset[ifc.PYRATE_X_STEP] = float(lookup[GAMMA_X_STEP][0])
    subset[ifc.PYRATE_DATUM] = "".join(lookup[GAMMA_DATUM])
    subset[ifc.PYRATE_INSAR_PROCESSOR] = GAMMA
    return subset


def _frequency_to_wavelength(freq):
    """
    Convert radar frequency to wavelength
    """
    return ifc.SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1, dem_hdr):
    """
    Combines metadata for master and slave epochs and DEM into a single
    dictionary for an interferogram.

    :param dict hdr0: Metadata for the master image
    :param dict hdr1: Metadata for the slave image
    :param dict dem_hdr: Metadata for the DEM

    :return: chdr: combined metadata
    :rtype: dict
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
            ifc.DATA_UNITS: RADIANS,
            ifc.PYRATE_INSAR_PROCESSOR: GAMMA}

    # set incidence angle to mean of master and slave
    inc_ang = hdr0[ifc.PYRATE_INCIDENCE_DEGREES]
    if np.isclose(inc_ang, hdr1[ifc.PYRATE_INCIDENCE_DEGREES], atol=1e-1):
        chdr[ifc.PYRATE_INCIDENCE_DEGREES] = (hdr0[ifc.PYRATE_INCIDENCE_DEGREES] + hdr1[
            ifc.PYRATE_INCIDENCE_DEGREES]) / 2
    else:
        msg = "Incidence angles differ by more than 1e-1"
        raise GammaException(msg)

    wavelen = hdr0[ifc.PYRATE_WAVELENGTH_METRES]
    if np.isclose(wavelen, hdr1[ifc.PYRATE_WAVELENGTH_METRES], atol=1e-6):
        chdr[ifc.PYRATE_WAVELENGTH_METRES] = wavelen
    else:
        args = (chdr[ifc.MASTER_DATE], chdr[ifc.SLAVE_DATE])
        msg = "Wavelength mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)
    # non-cropped, non-multilooked geotif process step information added
    chdr[ifc.DATA_TYPE] = ifc.ORIG

    chdr.update(dem_hdr)  # add geographic data
    return chdr


def manage_headers(dem_header_file, header_paths):
    """
    Manage and combine  header files for GAMMA interferograms, DEM and
    incidence files

    :param str dem_header_file: DEM header path
    :param list header_paths: List of master/slave header paths

    :return: combined_header: Combined metadata dictionary
    :rtype: dict
    """
    dem_header = parse_dem_header(dem_header_file)
    # find param files containing filename dates
    if len(header_paths) == 2:
        headers = [parse_epoch_header(hp) for hp in header_paths]
        combined_header = combine_headers(headers[0], headers[1], dem_header)
    else:
        # probably have DEM or incidence file
        combined_header = dem_header
        combined_header[ifc.DATA_TYPE] = ifc.DEM

    return combined_header


def get_header_paths(input_file, slc_file_list):
    """
    Function that matches input GAMMA file names with GAMMA header file names

    :param str input_file: input GAMMA image file.
    :param str slc_dir: GAMMA SLC header file directory
    :return: list of matching header files
    :rtype: list
    """
    f = Path(input_file)
    epochs = extract_epochs_from_filename(f.name)
    header_names = cf.parse_namelist(slc_file_list)
    matches = [hdr for hdr in header_names if any(e in hdr for e in epochs)]
    return matches


def gamma_header(ifg_file_path, params):
    """
    Function to obtain combined Gamma headers for image file
    
    Args:
        ifg_file_path: Path to interferogram file to find headers for.
        params: PyRate parameters dictionary.

    Returns:
        A combined header dictionary containing metadata from matching
        gamma headers and DEM header.   
    """
    dem_hdr_path = params[cf.DEM_HEADER_FILE]
    header_paths = get_header_paths(ifg_file_path, params[cf.HDR_FILE_LIST])
    combined_headers = manage_headers(dem_hdr_path, header_paths)
    if os.path.basename(ifg_file_path).split('.')[1] == \
            (params[cf.APS_INCIDENCE_EXT] or params[cf.APS_ELEVATION_EXT]):
        # TODO: implement incidence class here
        combined_headers['FILE_TYPE'] = 'Incidence'

    return combined_headers


class GammaException(Exception):
    """Gamma generic exception class"""
