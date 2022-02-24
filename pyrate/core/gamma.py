#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
from os.path import split
from pathlib import Path
from datetime import date, time, timedelta
import numpy as np

import pyrate.constants as C
from pyrate.configuration import ConfigException, parse_namelist
import pyrate.core.ifgconstants as ifc
from pyrate.constants import sixteen_digits_pattern, BASELINE_FILE_PATHS, BASE_FILE_DIR
from pyrate.core.shared import extract_epochs_from_filename, data_format
from pyrate.core.logger import pyratelogger as log
import struct


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
GAMMA_AZIMUTH = 'azimuth_angle'
GAMMA_RANGE_PIX = 'range_pixel_spacing'
GAMMA_RANGE_N = 'range_samples'
GAMMA_RANGE_LOOKS = 'range_looks'
GAMMA_AZIMUTH_PIX = 'azimuth_pixel_spacing'
GAMMA_AZIMUTH_N = 'azimuth_lines'
GAMMA_AZIMUTH_LOOKS = 'azimuth_looks'
GAMMA_PRF = 'prf'
GAMMA_NEAR_RANGE = 'near_range_slc'
GAMMA_SAR_EARTH = 'sar_to_earth_center'
GAMMA_SEMI_MAJOR_AXIS = 'earth_semi_major_axis'
GAMMA_SEMI_MINOR_AXIS = 'earth_semi_minor_axis'
GAMMA_INITIAL_BASELINE = 'initial_baseline(TCN)'
GAMMA_INITIAL_BASELINE_RATE = 'initial_baseline_rate'
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

    :param str path: Full path to GAMMA mli.par file. Note that the mli.par is required as input since the baseline calculations require the input values valid for the GAMMA multi-looked products and also the GAMMA lookup table gives radar coordinates for the multi-looked geometry.

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

    sat_azimuth, unit = lookup[GAMMA_AZIMUTH]
    if unit != "degrees":  # pragma: no cover
        msg = 'Unrecognised unit field for azimuth_angle: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_AZIMUTH_DEGREES] = float(sat_azimuth)

    range_pix, unit = lookup[GAMMA_RANGE_PIX]
    if unit != "m":  # pragma: no cover
        msg = 'Unrecognised unit field for range_pixel_spacing: %s'
        raise GammaException(msg % unit)
    subset[ifc.PYRATE_RANGE_PIX_METRES] = float(range_pix)

    range_n = lookup[GAMMA_RANGE_N] # number without a unit in .par file
    subset[ifc.PYRATE_RANGE_N] = int(range_n[0])

    range_looks = lookup[GAMMA_RANGE_LOOKS]  # number without a unit in .par file
    subset[ifc.PYRATE_RANGE_LOOKS] = int(range_looks[0])

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


def parse_baseline_header(path: str) -> dict:
    """
    Returns dictionary of Baseline metadata required for PyRate.
    Will read the Precise baseline estimate, if available,
    otherwise will read the Initial baseline estimate.

    :param path: Full path to GAMMA base.par file

    :return: bdict: Dictionary of baseline values
    """
    lookup = _parse_header(path)  # read file contents in to a dict

    # split the initial and precise baselines
    initial = lookup[GAMMA_INITIAL_BASELINE]
    initial_rate = lookup[GAMMA_INITIAL_BASELINE_RATE]
    precise = lookup[GAMMA_PRECISION_BASELINE]
    precise_rate = lookup[GAMMA_PRECISION_BASELINE_RATE]

    # read the initial baseline if all precise components are zero
    # (indicates that the precise baseline estimation was not ran in GAMMA workflow)
    if float(precise[0]) == 0.0 and float(precise[1]) == 0.0 and float(precise[2]) == 0.0:
        log.debug('Reading Initial GAMMA baseline values')
        baseline, baseline_rate = initial, initial_rate
    else:
        log.debug('Reading Precise GAMMA baseline values')
        baseline, baseline_rate = precise, precise_rate

    # Extract and return a dict of baseline values
    bdict = {}

    # baseline vector (along Track, aCross track, Normal to the track)
    bdict[ifc.PYRATE_BASELINE_T] = float(baseline[0])
    bdict[ifc.PYRATE_BASELINE_C] = float(baseline[1])
    bdict[ifc.PYRATE_BASELINE_N] = float(baseline[2])
    bdict[ifc.PYRATE_BASELINE_RATE_T] = float(baseline_rate[0])
    bdict[ifc.PYRATE_BASELINE_RATE_C] = float(baseline_rate[1])
    bdict[ifc.PYRATE_BASELINE_RATE_N] = float(baseline_rate[2])

    return bdict


def _frequency_to_wavelength(freq):
    """
    Convert radar frequency to wavelength
    """
    return ifc.SPEED_OF_LIGHT_METRES_PER_SECOND / freq


def combine_headers(hdr0, hdr1, dem_hdr, base_hdr=None):
    """
    Combines metadata for first and second image epochs, DEM and baselines
    into a single dictionary for an interferogram.

    :param dict hdr0: Metadata for the first image
    :param dict hdr1: Metadata for the second image
    :param dict dem_hdr: Metadata for the DEM
    :param dict base_hdr: Metadata for baselines (if available)

    :return: chdr: combined metadata
    :rtype: dict
    """
    if not all([isinstance(a, dict) for a in [hdr0, hdr1, dem_hdr]]):
        raise GammaException('Header args need to be dicts')

    if base_hdr and not isinstance(base_hdr, dict):
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
        msg = "Incidence angles differ by more than 0.1 degrees"
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
    if np.isclose(heading_ang, hdr1[ifc.PYRATE_HEADING_DEGREES], atol=5e-1):
        chdr[ifc.PYRATE_HEADING_DEGREES] = heading_ang
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Satellite heading angles differ by more than 0.5 degrees"
        raise GammaException(msg % args)

    azimuth_ang = hdr0[ifc.PYRATE_AZIMUTH_DEGREES]
    if np.isclose(azimuth_ang, hdr1[ifc.PYRATE_AZIMUTH_DEGREES], atol=1e-1):
        chdr[ifc.PYRATE_AZIMUTH_DEGREES] = azimuth_ang
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Satellite azimuth angles differ by more than 0.1 degrees"
        raise GammaException(msg % args)

    range_pix = hdr0[ifc.PYRATE_RANGE_PIX_METRES]
    if np.isclose(range_pix, hdr1[ifc.PYRATE_RANGE_PIX_METRES], atol=1e-1):
        chdr[ifc.PYRATE_RANGE_PIX_METRES] = range_pix
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Range pixel spacing differs by more than 0.001 metres"
        raise GammaException(msg % args)

    range_n = hdr0[ifc.PYRATE_RANGE_N]
    if range_n == hdr1[ifc.PYRATE_RANGE_N]:
        chdr[ifc.PYRATE_RANGE_N] = range_n
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Number of range pixels mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    range_looks = hdr0[ifc.PYRATE_RANGE_LOOKS]
    if range_looks == hdr1[ifc.PYRATE_RANGE_LOOKS]:
        chdr[ifc.PYRATE_RANGE_LOOKS] = range_looks
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Number of range looks mismatch, check both header files for %s & %s"
        raise GammaException(msg % args)

    azimuth_pix = hdr0[ifc.PYRATE_AZIMUTH_PIX_METRES]
    if np.isclose(azimuth_pix, hdr1[ifc.PYRATE_AZIMUTH_PIX_METRES], atol=1e-1):
        chdr[ifc.PYRATE_AZIMUTH_PIX_METRES] = azimuth_pix
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Azimuth pixel spacing differs by more than 0.001 metres"
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
        msg = "Near range differs by more than 1000 metres"
        raise GammaException(msg % args)

    sar_earth = hdr0[ifc.PYRATE_SAR_EARTH_METRES]
    if np.isclose(sar_earth, hdr1[ifc.PYRATE_SAR_EARTH_METRES], atol=1e3):
        chdr[ifc.PYRATE_SAR_EARTH_METRES] = sar_earth
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "SAR to Earth Center differs by more than 1000 metres"
        raise GammaException(msg % args)

    semi_major_axis = hdr0[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES]
    if np.isclose(semi_major_axis, hdr1[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES], atol=1e-4):
        chdr[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES] = semi_major_axis
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Earth semi major axis differs by more than 0.0001 metres"
        raise GammaException(msg % args)

    semi_minor_axis = hdr0[ifc.PYRATE_SEMI_MINOR_AXIS_METRES]
    if np.isclose(semi_minor_axis, hdr1[ifc.PYRATE_SEMI_MINOR_AXIS_METRES], atol=1e-4):
        chdr[ifc.PYRATE_SEMI_MINOR_AXIS_METRES] = semi_minor_axis
    else:
        args = (chdr[ifc.FIRST_DATE], chdr[ifc.SECOND_DATE])
        msg = "Earth semi minor axis differs by more than 0.0001 metres"
        raise GammaException(msg % args)

    # non-cropped, non-multilooked geotif process step information added
    chdr[ifc.DATA_TYPE] = ifc.ORIG

    chdr.update(dem_hdr)  # add geographic data

    if base_hdr:
        chdr.update(base_hdr)  # add baseline information

    return chdr


def manage_headers(dem_header_file, header_paths, baseline_paths=None):
    """
    Manage and combine header files for GAMMA interferograms, DEM and
    incidence files

    :param str dem_header_file: DEM header path
    :param list header_paths: List of first/second image header paths

    :return: combined_header: Combined metadata dictionary
    :rtype: dict
    """
    dem_header = parse_dem_header(dem_header_file)
    # find param files containing filename dates
    if len(header_paths) == 2:
        hdrs = [parse_epoch_header(hp) for hp in header_paths]
        if baseline_paths is not None:
            baseline_header = parse_baseline_header(baseline_paths)
            combined_header = combine_headers(hdrs[0], hdrs[1], dem_header, baseline_header)
        else:
            combined_header = combine_headers(hdrs[0], hdrs[1], dem_header)
    elif len(header_paths) > 2:
        msg = f'There are too many parameter files for one interferogram; there should only be two. {len(header_paths)} parameter files have been given: {header_paths}.'
        raise GammaException(msg)
    else:
        # probably have DEM or incidence file
        combined_header = dem_header
        combined_header[ifc.DATA_TYPE] = ifc.DEM

    return combined_header


def get_header_paths(input_file, slc_file_list):
    """
    Function that matches input GAMMA file names with GAMMA header file names

    :param str input_file: input GAMMA image file.
    :param slc_file_list: file listing the pool of available header files (GAMMA: slc.par, ROI_PAC: .rsc)
    :return: list of matching header files
    :rtype: list
    """
    f = Path(input_file)
    epochs = extract_epochs_from_filename(f.name)
    header_names = parse_namelist(slc_file_list)
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
    dem_hdr_path = params[C.DEM_HEADER_FILE]
    header_paths = get_header_paths(ifg_file_path, params[C.HDR_FILE_LIST])
    if len(header_paths) == 2 and params[C.BASE_FILE_LIST] is not None:
        baseline_path = baseline_paths_for(ifg_file_path, params)
    else:
        baseline_path = None  # don't read baseline files for DEM

    combined_headers = manage_headers(dem_hdr_path, header_paths, baseline_path)

    if os.path.basename(ifg_file_path).split('.')[1] == \
            (params[C.APS_INCIDENCE_EXT] or params[C.APS_ELEVATION_EXT]):
        # TODO: implement incidence class here
        combined_headers['FILE_TYPE'] = 'Incidence'

    return combined_headers


def read_lookup_table(head, data_path, xlooks, ylooks, xmin, xmax, ymin, ymax):
    # pylint: disable = too - many - statements
    """
    Creates a copy of input lookup table file in a numpy array and applies the ifg ML factors

    :param IFG object head: first IFG in the list to read metadata
    :param str data_path: Input file
    :param int xlooks: multi-looking factor in x
    :param int ylooks: multi-looking factor in y
    :param int xmin: start pixel of cropped extent in x
    :param int xmax: end pixel of cropped extent in x
    :param int ymin: start pixel of cropped extent in y
    :param int ymax: end pixel of cropped extent in y

    :return: np-array lt_data_az: azimuth (i.e. row) of radar-coded MLI
    :return: np-array lt_data_rg: range (i.e. column) of radar-coded MLI
    """
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-locals

    # read relevant metadata parameters
    nrows_lt = int(head.meta_data[ifc.PYRATE_NROWS]) # number of rows of original geotiff files
    ncols_lt = int(head.meta_data[ifc.PYRATE_NCOLS]) # number of columns of original geotiff files
    nrows = head.nrows # number of rows in multi-looked and cropped data sets
    ncols = head.ncols # number of columns in multi-looked and cropped data sets

    ifg_proc = head.meta_data[ifc.PYRATE_INSAR_PROCESSOR]
    # get dimensions of lookup table file
    bytes_per_col, fmtstr = data_format(ifg_proc, True, ncols_lt*2) # float complex data set containing value tupels

    # check if lookup table has the correct size
    small_size = _check_raw_data(bytes_per_col * 2, data_path, ncols_lt, nrows_lt)
    # todo: delete the following if condition once a suitable test data set has been included
    if small_size: # this is a test data set without a corresponding lt-file
        lt_data_az = np.empty((nrows, ncols)) * np.nan # nan array with size of input data set
        lt_data_rg = np.empty((nrows, ncols)) * np.nan # nan array with size of input data set

    else: # this is a real data set with an lt-file of correct size
        row_bytes = ncols_lt * 2 * bytes_per_col
        lt_data_az = np.empty((0, ncols))  # empty array with correct number of columns
        lt_data_rg = np.empty((0, ncols))  # empty array with correct number of column
        # for indexing: lookup table file contains value pairs (i.e. range, azimuth)
        # value pair 0 would be index 0 and 1, value pair 1 would be index 2 and 3, and so on
        # example: for a multi-looking factor of 10 we want value pair 4, 14, 24, ...
        # this would be index 8 and 9, index 28 and 29, 48 and 49, ...
        # start column needs to be added in case cropping is applied
        if (xlooks % 2) == 0:  # for even ml factors
            idx_start = xmin + int(xlooks / 2) - 1
        else: # for odd ml factors
            idx_start = xmin + int((xlooks - 1) / 2)
        # indices of range info in lookup table for the cropped and multi-looked data set
        idx_rg = np.arange(2 * idx_start, xmax * 2, 2 * xlooks)  # first value
        idx_az = np.arange(2 * idx_start + 1, xmax * 2, 2 * xlooks)  # second value
        # set up row idx, e.g. for ml=10 (without cropping): 4, 14, 24, ...
        if (ylooks % 2) == 0: # for even ml factors
            idx_start = ymin + int(ylooks / 2) - 1
        else: # for odd ml factors
            idx_start = ymin + int((ylooks - 1) / 2)
        row_idx = np.arange(idx_start, ymax, ylooks)

        # read the binary lookup table file and save the range/azimuth value pair for each position in the cropped and
        # multi-looked data set
        log.debug(f"Reading lookup table file {data_path}")
        with open(data_path, 'rb') as f:
            for y in range(nrows_lt): # loop through all lines in file
                # this could potentially be made quicker by skipping unwanted bytes in the f.read command?
                data = struct.unpack(fmtstr, f.read(row_bytes))
                # but only read data from lines in row index:
                if y in row_idx:
                    row_data = np.array(data)
                    row_data_ml_az = row_data[idx_az] # azimuth for PyRate
                    row_data_ml_rg = row_data[idx_rg] # range for PyRate
                    lt_data_az = np.append(lt_data_az, [row_data_ml_az], axis=0)
                    lt_data_rg = np.append(lt_data_rg, [row_data_ml_rg], axis=0)

    return lt_data_az, lt_data_rg


def _check_raw_data(bytes_per_col, data_path, ncols, nrows):
    """
    Convenience function to check the file size is as expected
    """
    size = ncols * nrows * bytes_per_col
    act_size = os.stat(data_path).st_size
    if act_size != size:
        msg = '%s should have size %s, not %s. Is the correct file being used?'
        if size < 28000:
            # test data set doesn't currently fit the lookup table size, stop further calculation
            # todo: delete this if statement once a new test data set has been introduced
            return True
        else:
            raise GammaException(msg % (data_path, size, act_size))


class GammaException(Exception):
    """Gamma generic exception class"""


def baseline_paths_for(path: str, params: dict) -> str:
    """
    Returns path to baseline file for given interferogram. Pattern matches
    based on epoch in filename.

    Example:
        '20151025-20160501_base.par'
        Date pair is the epoch.

    Args:
        path: Path to intergerogram to find baseline file for.
        params: Parameter dictionary.
        tif: Find converted tif if True (_cc.tif), else find .cc file.

    Returns:
        Path to baseline file.
    """

    _, filename = split(path)
    try:
        epoch = re.search(sixteen_digits_pattern, filename).group(0)
    except: # catch cases where filename does not have two epochs, e.g. DEM file
        return None

    base_file_paths = [f.unwrapped_path for f in params[BASELINE_FILE_PATHS] if epoch in f.unwrapped_path]

    if len(base_file_paths) > 1:
        raise ConfigException(f"'{BASE_FILE_DIR}': found more than one baseline "
                              f"file for '{path}'. There must be only one "
                              f"baseline file per interferogram. Found {base_file_paths}.")
    return base_file_paths[0]
