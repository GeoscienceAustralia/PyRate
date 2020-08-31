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
GAMMA_HEADING = 'heading'
GAMMA_RANGE_PIX = 'range_pixel_spacing'
GAMMA_RANGE_N = 'range_samples'
GAMMA_AZIMUTH_PIX = 'azimuth_pixel_spacing'
GAMMA_AZIMUTH_N = 'azimuth_lines'
GAMMA_AZIMUTH_LOOKS = 'azimuth_looks'
GAMMA_PRF = 'prf'
GAMMA_NEAR_RANGE = 'near_range_slc'
GAMMA_SAR_EARTH = 'sar_to_earth_center'
GAMMA_SEMI_MAJOR_AXIS = 'earth_semi_major_axis'
GAMMA_SEMI_MINOR_AXIS = 'earth_semi_minor_axis'
# to do add option to use initial baseline if precision baseline was not calculated in GAMMA
GAMMA_PRECISION_BASELINE = 'precision_baseline(TCN)'
GAMMA_PRECISION_BASELINE_RATE = 'precision_baseline_rate'
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

    sat_heading, unit = lookup[GAMMA_HEADING]
    if unit != "degrees":  # pragma: no cover
        msg = 'Unrecognised unit field for heading: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_HEADING_DEGREES] = float(sat_heading)

    range_pix, unit = lookup[GAMMA_RANGE_PIX]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for range_pixel_spacing: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_RANGE_PIX_METRES] = float(range_pix)

    range_n = lookup[GAMMA_RANGE_N] # number without a unit in .par file
    subset[ifc.PYRATE_RANGE_N] = int(range_n[0])

    azimuth_pix, unit = lookup[GAMMA_AZIMUTH_PIX]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for azimuth_pixel_spacing: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_AZIMUTH_PIX_METRES] = float(azimuth_pix)

    azimuth_n = lookup[GAMMA_AZIMUTH_N] # number without a unit in .par file
    subset[ifc.PYRATE_AZIMUTH_N] = int(azimuth_n[0])

    azimuth_looks = lookup[GAMMA_AZIMUTH_LOOKS]  # number without a unit in .par file
    subset[ifc.PYRATE_AZIMUTH_LOOKS] = int(azimuth_looks[0])

    pulse_rep_freq, unit = lookup[GAMMA_PRF]
    if unit != "Hz":  # pragma: no cover
        msg = 'Unrecognised unit field for prf: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_PRF_HERTZ] = float(pulse_rep_freq)

    near_range, unit = lookup[GAMMA_NEAR_RANGE]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for near_range_slc: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_NEAR_RANGE_METRES] = float(near_range)

    sar_to_earth, unit = lookup[GAMMA_SAR_EARTH]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for sar_to_earth_center: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_SAR_EARTH_METRES] = float(sar_to_earth)

    semi_major_axis, unit = lookup[GAMMA_SEMI_MAJOR_AXIS]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for earth_semi_major_axis: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES] = float(semi_major_axis)

    semi_minor_axis, unit = lookup[GAMMA_SEMI_MINOR_AXIS]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for earth_semi_minor_axis: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_SEMI_MINOR_AXIS_METRES] = float(semi_minor_axis)

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

    subset[ifc.FIRST_DATE] = date(year, month, day)
    subset[ifc.FIRST_TIME] = time(hour, min, sec)

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


def parse_baseline_header(path):
    """
    Returns dictionary of DEM metadata required for PyRate

    :param str path: `Full path to Gamma *dem.par file`

    :return: subset: subset of full metadata
    :rtype: dict
    """
    lookup = _parse_header(path)

    subset = {}
    # baseline vector (along Track, aCross track, Normal to the track)
    baseline_tcn = lookup[GAMMA_PRECISION_BASELINE]
    subset[ifc.PYRATE_BASELINE_T] = float(baseline_tcn[0])
    subset[ifc.PYRATE_BASELINE_C] = float(baseline_tcn[1])
    subset[ifc.PYRATE_BASELINE_N] = float(baseline_tcn[2])
    baseline_rate_tcn = lookup[GAMMA_PRECISION_BASELINE_RATE]
    subset[ifc.PYRATE_BASELINE_RATE_T] = float(baseline_rate_tcn[0])
    subset[ifc.PYRATE_BASELINE_RATE_C] = float(baseline_rate_tcn[1])
    subset[ifc.PYRATE_BASELINE_RATE_N] = float(baseline_rate_tcn[2])

    return subset


def _frequency_to_wavelength(freq):
    """
    Convert radar frequency to wavelength
    """
    return ifc.SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1, dem_hdr, bas_hdr):
    """
    Combines metadata for first and second image epochs and DEM into a
    single dictionary for an interferogram.

    :param dict hdr0: Metadata for the first image
    :param dict hdr1: Metadata for the second image
    :param dict dem_hdr: Metadata for the DEM

    :return: chdr: combined metadata
    :rtype: dict
    """
    if not all([isinstance(a, dict) for a in [hdr0, hdr1, dem_hdr, bas_hdr]]):
        raise GammaException('Header args need to be dicts')

    date0, date1 = hdr0[ifc.FIRST_DATE], hdr1[ifc.FIRST_DATE]
    if date0 == date1:
        raise GammaException("Can't combine headers for the same day")
    elif date1 < date0:
        raise GammaException("Wrong date order")

    chdr = {ifc.PYRATE_TIME_SPAN: (date1 - date0).days / ifc.DAYS_PER_YEAR,
            ifc.FIRST_DATE: date0,
            ifc.FIRST_TIME: hdr0[ifc.FIRST_TIME],
            ifc.SECOND_DATE: date1,
            ifc.SECOND_TIME: hdr1[ifc.FIRST_TIME],
            ifc.DATA_UNITS: RADIANS,
            ifc.PYRATE_INSAR_PROCESSOR: GAMMA}

    # set incidence angle to mean of first amd second image values
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
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Wavelength mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    # use parameter of first image (as done by GAMMA during interferometric processing)
    heading_ang = hdr0[ifc.PYRATE_HEADING_DEGREES]
    if np.isclose(heading_ang, hdr1[ifc.PYRATE_HEADING_DEGREES], atol=1e-1):
        chdr[ifc.PYRATE_HEADING_DEGREES] = heading_ang
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Satellite heading angles differ by more than 1e-1"
        raise GammaException(msg % args)

    range_pix = hdr0[ifc.PYRATE_RANGE_PIX_METRES]
    if np.isclose(range_pix, hdr1[ifc.PYRATE_RANGE_PIX_METRES], atol=1e-1):
        chdr[ifc.PYRATE_RANGE_PIX_METRES] = range_pix
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Range pixel spacing differs by more than 1e-3"
        raise GammaException(msg % args)

    range_n = hdr0[ifc.PYRATE_RANGE_N]
    if range_n == hdr1[ifc.PYRATE_RANGE_N]:
        chdr[ifc.PYRATE_RANGE_N] = range_n
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Number of range pixels mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    azimuth_pix = hdr0[ifc.PYRATE_AZIMUTH_PIX_METRES]
    if np.isclose(azimuth_pix, hdr1[ifc.PYRATE_AZIMUTH_PIX_METRES], atol=1e-1):
        chdr[ifc.PYRATE_AZIMUTH_PIX_METRES] = azimuth_pix
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Azimuth pixel spacing differs by more than 1e-3"
        raise GammaException(msg % args)

    azimuth_n = hdr0[ifc.PYRATE_AZIMUTH_N]
    if azimuth_n == hdr1[ifc.PYRATE_AZIMUTH_N]:
        chdr[ifc.PYRATE_AZIMUTH_N] = azimuth_n
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Number of azimuth pixels mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    azimuth_looks = hdr0[ifc.PYRATE_AZIMUTH_LOOKS]
    if azimuth_looks == hdr1[ifc.PYRATE_AZIMUTH_LOOKS]:
        chdr[ifc.PYRATE_AZIMUTH_LOOKS] = azimuth_looks
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Number of azimuth looks mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    prf_hertz = hdr0[ifc.PYRATE_PRF_HERTZ]
    if np.isclose(prf_hertz, hdr1[ifc.PYRATE_PRF_HERTZ], atol=1e-6):
        chdr[ifc.PYRATE_PRF_HERTZ] = prf_hertz
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Pulse repetition frequency mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    near_range = hdr0[ifc.PYRATE_NEAR_RANGE_METRES]
    if np.isclose(near_range, hdr1[ifc.PYRATE_NEAR_RANGE_METRES], atol=1e3):
        chdr[ifc.PYRATE_NEAR_RANGE_METRES] = near_range
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Near range differs by more than 1e3"
        raise GammaException(msg % args)

    sar_earth = hdr0[ifc.PYRATE_SAR_EARTH_METRES]
    if np.isclose(sar_earth, hdr1[ifc.PYRATE_SAR_EARTH_METRES], atol=1e3):
        chdr[ifc.PYRATE_SAR_EARTH_METRES] = sar_earth
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "SAR to Earth Center differs by more than 1e3"
        raise GammaException(msg % args)

    semi_major_axis = hdr0[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES]
    if np.isclose(semi_major_axis, hdr1[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES], atol=1e-4):
        chdr[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES] = semi_major_axis
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Earth semi major axis differs by more than 1e-4"
        raise GammaException(msg % args)

    semi_minor_axis = hdr0[ifc.PYRATE_SEMI_MINOR_AXIS_METRES]
    if np.isclose(semi_minor_axis, hdr1[ifc.PYRATE_SEMI_MINOR_AXIS_METRES], atol=1e-4):
        chdr[ifc.PYRATE_SEMI_MINOR_AXIS_METRES] = semi_minor_axis
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Earth semi minor axis differs by more than 1e-4"
        raise GammaException(msg % args)

    # non-cropped, non-multilooked geotif process step information added
    chdr[ifc.DATA_TYPE] = ifc.ORIG

    chdr.update(dem_hdr)  # add geographic data

    chdr.update(bas_hdr)  # add baseline information

    return chdr


def manage_headers(dem_header_file, header_paths, baseline_paths):
    """
    Manage and combine  header files for GAMMA interferograms, DEM and
    incidence files

    :param str dem_header_file: DEM header path
    :param list header_paths: List of first/second image header paths

    :return: combined_header: Combined metadata dictionary
    :rtype: dict
    """
    dem_header = parse_dem_header(dem_header_file)
    # find param files containing filename dates
    if len(header_paths) == 2:
        headers = [parse_epoch_header(hp) for hp in header_paths]
        baseline_header = parse_baseline_header(baseline_paths)
        combined_header = combine_headers(headers[0], headers[1], dem_header, baseline_header)
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
    print(ifg_file_path)
    if len(header_paths) == 2:
        baseline_path = cf.baseline_paths_for(ifg_file_path, params)
    else:
        baseline_path = ''
    combined_headers = manage_headers(dem_hdr_path, header_paths, baseline_path)
    if os.path.basename(ifg_file_path).split('.')[1] == \
            (params[cf.APS_INCIDENCE_EXT] or params[cf.APS_ELEVATION_EXT]):
        # TODO: implement incidence class here
        combined_headers['FILE_TYPE'] = 'Incidence'

    return combined_headers


class GammaException(Exception):
    """Gamma generic exception class"""
