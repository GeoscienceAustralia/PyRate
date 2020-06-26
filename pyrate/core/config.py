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
This Python module contains utilities to parse PyRate configuration
files. It also includes numerous general constants relating to options
in configuration files.
"""
# coding: utf-8
# pylint: disable=invalid-name
# pylint: disable=W1203
# pylint: disable=too-many-locals
# pylint: disable=trailing-whitespace
from typing import List, Tuple, Dict, Optional
import os
from os.path import splitext, split
import re
from pathlib import Path
from osgeo import gdal

from pyrate.core.ifgconstants import YEARS_PER_DAY
from pyrate.constants import CONV2TIF, PREPIFG, PROCESS, MERGE
from pyrate.constants import SIXTEEN_DIGIT_EPOCH_PAIR, TWELVE_DIGIT_EPOCH_PAIR, EIGHT_DIGIT_EPOCH, \
    sixteen_digits_pattern
from pyrate.core.logger import pyratelogger as _logger

# general constants
MINIMUM_NUMBER_EPOCHS = 3
NO_MULTILOOKING = 1
ROIPAC = 0
GAMMA = 1
LOG_LEVEL = 'INFO'

# constants for lookups
#: STR; Name of input interferogram list file
IFG_FILE_LIST = 'ifgfilelist'
#: (0/1/2); The interferogram processor used (0==ROIPAC, 1==GAMMA, 2: GEOTIF)
PROCESSOR = 'processor'
#: STR; Name of directory containing input interferograms.
OBS_DIR = 'obsdir'
#: STR; Name of directory for saving output products
OUT_DIR = 'outdir'
#: STR; Name of Digital Elevation Model file
DEM_FILE = 'demfile'
#: STR; Name of the DEM header file
DEM_HEADER_FILE = 'demHeaderFile'
#: STR; Name of directory containing GAMMA SLC header files
SLC_DIR = 'slcFileDir'
#: STR; Name of the file list containing the pool of available header files
HDR_FILE_LIST = 'hdrfilelist'


INTERFEROGRAM_FILES = 'interferogram_files'
HEADER_FILE_PATHS = 'header_file_paths'
COHERENCE_FILE_PATHS = 'coherence_file_paths'
DEM_FILE_PATH = 'dem_file'

# STR; The projection of the input interferograms.
# TODO: only used in tests; deprecate?
INPUT_IFG_PROJECTION = 'projection'
#: FLOAT; The no data value in the interferogram files.
NO_DATA_VALUE = 'noDataValue'
#: FLOAT; No data averaging threshold for prepifg
NO_DATA_AVERAGING_THRESHOLD = 'noDataAveragingThreshold'
# BOOL (1/2/3); Re-project data from Line of sight, 1 = vertical, 2 = horizontal, 3 = no conversion
#REPROJECTION = 'prjflag' # NOT CURRENTLY USED
#: BOOL (0/1): Convert no data values to Nan
NAN_CONVERSION = 'nan_conversion'

# Prepifg parameters
#: BOOL (1/2/3/4); Method for cropping interferograms, 1 = minimum overlapping area (intersection), 2 = maximum area (union), 3 = customised area, 4 = all ifgs already same size
IFG_CROP_OPT = 'ifgcropopt'
#: INT; Multi look factor for interferogram preparation in x dimension
IFG_LKSX = 'ifglksx'
#: INT; Multi look factor for interferogram preparation in y dimension
IFG_LKSY = 'ifglksy'
#: FLOAT; Minimum longitude for cropping with method 3
IFG_XFIRST = 'ifgxfirst'
#: FLOAT; Maximum longitude for cropping with method 3
IFG_XLAST = 'ifgxlast'
#: FLOAT; Minimum latitude for cropping with method 3
IFG_YFIRST = 'ifgyfirst'
#: FLOAT; Maximum latitude for cropping with method 3
IFG_YLAST = 'ifgylast'

# reference pixel parameters
#: INT; Longitude (decimal degrees) of reference pixel, or if left blank a search will be performed
REFX = 'refx'
#: INT; Latitude (decimal degrees) of reference pixel, or if left blank a search will be performed
REFY = 'refy'
#: INT; Number of reference pixel grid search nodes in x dimension
REFNX = "refnx"
#: INT; Number of reference pixel grid search nodes in y dimension
REFNY = "refny"
#: INT; Dimension of reference pixel search window (in number of pixels)
REF_CHIP_SIZE = 'refchipsize'
#: FLOAT; Minimum fraction of observations required in search window for pixel to be a viable reference pixel
REF_MIN_FRAC = 'refminfrac'
#: BOOL (1/2); Reference phase estimation method (1: median of the whole interferogram, 2: median within the window surrounding the reference pixel)
REF_EST_METHOD = 'refest'

# coherence masking parameters
#: BOOL (0/1); Perform coherence masking (1: yes, 0: no)
COH_MASK = 'cohmask'
#: FLOAT; Coherence threshold for masking
COH_THRESH = 'cohthresh'
#: STR; Directory containing coherence files; defaults to OBS_DIR if not provided
COH_FILE_DIR = 'cohfiledir'
#: STR; Name of the file list containing the pool of available coherence files
COH_FILE_LIST = 'cohfilelist'

#atmospheric error correction parameters NOT CURRENTLY USED
APS_CORRECTION = 'apscorrect'
APS_METHOD = 'apsmethod'
APS_INCIDENCE_MAP = 'incidencemap'
APS_INCIDENCE_EXT = 'APS_INCIDENCE_EXT'
APS_ELEVATION_MAP = 'elevationmap'
APS_ELEVATION_EXT = 'APS_ELEVATION_EXT'

# orbital error correction/parameters
#: BOOL (1/0); Perform orbital error correction (1: yes, 0: no)
ORBITAL_FIT = 'orbfit'
#: BOOL (1/2); Method for orbital error correction (1: independent, 2: network)
ORBITAL_FIT_METHOD = 'orbfitmethod'
#: BOOL (1/2/3) Polynomial order of orbital error model (1: planar in x and y - 2 parameter model, 2: quadratic in x and y - 5 parameter model, 3: quadratic in x and cubic in y - part-cubic 6 parameter model)
ORBITAL_FIT_DEGREE = 'orbfitdegrees'
#: INT; Multi look factor for orbital error calculation in x dimension
ORBITAL_FIT_LOOKS_X = 'orbfitlksx'
#: INT; Multi look factor for orbital error calculation in y dimension
ORBITAL_FIT_LOOKS_Y = 'orbfitlksy'

# Stacking parameters
#: FLOAT; Threshold ratio between 'model minus observation' residuals and a-priori observation standard deviations for stacking estimate acceptance (otherwise remove furthest outlier and re-iterate)
LR_NSIG = 'nsig'
#: INT; Number of required observations per pixel for stacking to occur
LR_PTHRESH = 'pthr'
#: FLOAT; Maximum allowable standard error for pixels in stacking
LR_MAXSIG = 'maxsig'

# atmospheric delay errors fitting parameters NOT CURRENTLY USED
# atmfitmethod = 1: interferogram by interferogram; atmfitmethod = 2, epoch by epoch
#ATM_FIT = 'atmfit'
#ATM_FIT_METHOD = 'atmfitmethod'

# APS correction parameters
#: BOOL (0/1) Perform APS correction (1: yes, 0: no)
APSEST = 'apsest'
# temporal low-pass filter parameters
#: INT (1/2/3); Method for temporal filtering (1: Gaussian, 2: Triangular, 3: Mean filter)
TLPF_METHOD = 'tlpfmethod'
#: FLOAT; Cutoff time for gaussian filter in years;
TLPF_CUTOFF = 'tlpfcutoff'
#: INT; Number of required input observations per pixel for temporal filtering
TLPF_PTHR = 'tlpfpthr'
# spatially correlated noise low-pass filter parameters
#: INT (1/2); Method for spatial filtering(1: butterworth; 2: gaussian)
SLPF_METHOD = 'slpfmethod'
#: FLOAT; Cutoff  value for both butterworth and gaussian filters in km
SLPF_CUTOFF = 'slpfcutoff'
#: INT; Order of butterworth filter (default 1)
SLPF_ORDER = 'slpforder'
#: INT (1/0); Do spatial interpolation at NaN locations (1 for interpolation, 0 for zero fill)
SLPF_NANFILL = 'slpnanfill'
#: #: STR; Method for spatial interpolation (one of: linear, nearest, cubic), only used when slpnanfill=1
SLPF_NANFILL_METHOD = 'slpnanfill_method'

# Time series parameters
#: BOOL (1/0); Perform time series calculation (1: yes, 0: no)
TIME_SERIES_CAL = 'tscal'
#: INT (1/2); Method for time series inversion (1: Laplacian Smoothing; 2: SVD)
TIME_SERIES_METHOD = 'tsmethod'
#: INT; Number of required input observations per pixel for time series inversion
TIME_SERIES_PTHRESH = 'ts_pthr'
#: INT (1/2); Order of Laplacian smoothing operator, first or second order
TIME_SERIES_SM_ORDER = 'smorder'
#: FLOAT; Laplacian smoothing factor (values used is 10**smfactor)
TIME_SERIES_SM_FACTOR = 'smfactor'
# tsinterp is automatically assigned in the code; not needed in conf file
#TIME_SERIES_INTERP = 'tsinterp'

#: BOOL (0/1/2); Use parallelisation/Multi-threading (0: in serial, 1: in parallel by rows, 2: in parallel by pixel)
PARALLEL = 'parallel'
#: INT; Number of processes for multi-threading
PROCESSES = 'processes'
LARGE_TIFS = 'largetifs'
# Orbital error correction constants for conversion to readable strings
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2
PLANAR = 1
QUADRATIC = 2
PART_CUBIC = 3

# Orbital error name look up for logging
ORB_METHOD_NAMES = {INDEPENDENT_METHOD: 'INDEPENDENT', 
                    NETWORK_METHOD: 'NETWORK'}
ORB_DEGREE_NAMES = {PLANAR: 'PLANAR', 
                    QUADRATIC: 'QUADRATIC', 
                    PART_CUBIC: 'PART CUBIC'}

# dir for temp files
TMPDIR = 'tmpdir'

# Lookup to help convert args to correct type/defaults
# format is	key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = {
#    REPROJECTION : (int, 3), # Default no conversion, CONVERSION NOT IMPLEMENTED
    IFG_CROP_OPT : (int, 1), # default to area 'intersection' option
    IFG_LKSX : (int, NO_MULTILOOKING),
    IFG_LKSY : (int, NO_MULTILOOKING),
    IFG_XFIRST : (float, None),
    IFG_XLAST : (float, None),
    IFG_YFIRST : (float, None),
    IFG_YLAST : (float, None),
    NO_DATA_VALUE: (float, 0.0),

    COH_MASK: (int, 0),
    COH_THRESH: (float, 0.1),

    REFX: (float, -1),
    REFY: (float, -1),
    REFNX: (int, 10),
    REFNY: (int, 10),
    REF_CHIP_SIZE: (int, 21),
    REF_MIN_FRAC: (float, 0.5),
    REF_EST_METHOD: (int, 1),  # default to average of whole image

    ORBITAL_FIT: (int, 0),
    ORBITAL_FIT_METHOD: (int, NETWORK_METHOD),
    ORBITAL_FIT_DEGREE: (int, PLANAR),
    ORBITAL_FIT_LOOKS_X: (int, 10),
    ORBITAL_FIT_LOOKS_Y: (int, 10),

    LR_NSIG: (int, 2),
    # pixel thresh based on nepochs? not every project may have 20 epochs
    LR_PTHRESH: (int, 3),
    LR_MAXSIG: (int, 10),

    #ATM_FIT: (int, 0), NOT CURRENTLY USED
    #ATM_FIT_METHOD: (int, 2),

    APSEST: (int, 0),
    TLPF_METHOD: (int, 1),
    TLPF_CUTOFF: (float, 1.0),
    TLPF_PTHR: (int, 1),

    SLPF_METHOD: (int, 1),
    SLPF_CUTOFF: (float, 1.0),
    SLPF_ORDER: (int, 1),
    SLPF_NANFILL: (int, 0),

    TIME_SERIES_CAL: (int, 0),
    # pixel thresh based on nepochs? not every project may have 20 epochs
    TIME_SERIES_PTHRESH: (int, 3),
    TIME_SERIES_SM_FACTOR: (float, -1.0),
    TIME_SERIES_SM_ORDER: (int, None),
    TIME_SERIES_METHOD: (int, 2),  # Default to SVD method

    PARALLEL: (int, 0),
    PROCESSES: (int, 8),
    PROCESSOR: (int, None),
    NAN_CONVERSION: (int, 0),
    NO_DATA_AVERAGING_THRESHOLD: (float, 0.0),
    }

PATHS = [
    OBS_DIR,
    IFG_FILE_LIST,
    DEM_FILE,
    DEM_HEADER_FILE,
    OUT_DIR,
    SLC_DIR,
    HDR_FILE_LIST,
    COH_FILE_DIR,
    COH_FILE_LIST,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
]

DEFAULT_TO_OBS_DIR = [SLC_DIR, COH_FILE_DIR]

INT_KEYS = [APS_CORRECTION, APS_METHOD]


def get_config_params(path: str) -> Dict:
    """
    Reads the parameters file provided by the user and converts it into
    a dictionary.

    Args:
        path: Absolute path to the parameters file.
        validate: Validate the parameters if True, otherwise skip validation.
        step: The current step of the PyRate workflow.

    Returns:
       A dictionary of parameters.
    """
    txt = ''
    with open(path, 'r') as inputFile:
        for line in inputFile:
            if any(x in line for x in PATHS):
                pos = line.find('~')
                if pos != -1:
                    # create expanded line
                    line = line[:pos] + os.environ['HOME'] + line[(pos+1):]
            txt += line
    params = _parse_conf_file(txt)
    params[TMPDIR] = os.path.join(os.path.abspath(params[OUT_DIR]), 'tmpdir')

    return params


def _parse_conf_file(content) -> Dict:
    """
    Converts the parameters from their text form into a dictionary.

    Args:
        content: Parameters as text.

    Returns:
        A dictionary of parameters.
    """
    def _is_valid(line):
        """
        Check if line is not empty or has % or #
        """
        return line != "" and line[0] not in "%#"

    lines = [ln.split() for ln in content.split('\n') if _is_valid(ln)]

    # convert "field:   value" lines to [field, value]
    kvpair = [(e[0].rstrip(":"), e[1]) for e in lines if len(e) == 2] + \
             [(e[0].rstrip(":"), None) for e in lines if len(e) == 1]
    parameters = dict(kvpair)
    for p in PATHS:
        if p not in parameters:
            parameters[p] = None

    for p in INT_KEYS:
        if p not in parameters:
            parameters[p] = '0'  # insert dummies

    parameters = _handle_extra_parameters(parameters)

    if not parameters:
        raise ConfigException('Cannot parse any parameters from config file')

    return _parse_pars(parameters)


def _handle_extra_parameters(params):
    """
    Function to check if requirements for weather model correction are given.
    """
    params[APS_INCIDENCE_EXT] = None
    params[APS_ELEVATION_EXT] = None

    if params[APS_INCIDENCE_MAP] is not None:
        params[APS_INCIDENCE_EXT] = \
            os.path.basename(params[APS_INCIDENCE_MAP]).split('.')[-1]
        params[APS_ELEVATION_MAP] = None
        params[APS_ELEVATION_EXT] = None
        return params

    # define APS_ELEVATON_EXT for gamma prepifg
    if params[APS_ELEVATION_MAP] is not None:
        params[APS_ELEVATION_EXT] = os.path.basename(
            params[APS_ELEVATION_MAP]).split('.')[-1]

    return params


def _parse_pars(pars) -> Dict:
    """
    Takes dictionary of parameters, converting values to required type
    and providing defaults for missing values.

    Args:
        pars: Parameters dictionary.

    Returns:
        Dictionary of converted (and optionally validated) parameters.
    """
    # Fallback to default for missing values and perform conversion.
    for k in PARAM_CONVERSION:
        if pars.get(k) is None:
            pars[k] = PARAM_CONVERSION[k][1]
            # _logger.warning(f"No value found for parameter '{k}'. Using "f"default value {pars[k]}.")
        else:
            conversion_func = PARAM_CONVERSION[k][0]
            if conversion_func:
                try:
                    pars[k] = conversion_func(pars[k])
                except ValueError as e:
                    _logger.error(
                        f"Unable to convert '{k}': {pars[k]} to " f"expected type {conversion_func.__name__}.")
                    raise e

    # Fallback to default for missing paths.
    for p in DEFAULT_TO_OBS_DIR:
        if pars.get(p) is None:
            pars[p] = pars[OBS_DIR]

    return pars


# CONFIG UTILS - TO BE MOVED?
def parse_namelist(nml):
    """
    Parses name list file into array of paths

    :param str nml: interferogram file list

    :return: list of interferogram file names
    :rtype: list
    """
    with open(nml) as f_in:
        lines = [line.rstrip() for line in f_in]
    return filter(None, lines)


def write_config_file(params, output_conf_file):
    """
    Takes a param object and writes the config file. Reverse of get_conf_params.

    :param dict params: parameter dictionary
    :param str output_conf_file: output file name

    :return: config file
    :rtype: list
    """
    with open(output_conf_file, 'w') as f:
        for k, v in params.items():
            if v is not None:
                f.write(''.join([k, ':\t', str(v), '\n']))
            else:
                f.write(''.join([k, ':\t', '', '\n']))


def transform_params(params):
    """
    Returns subset of all parameters for cropping and multilooking.

    :param dict params: Parameter dictionary

    :return: xlooks, ylooks, crop
    :rtype: int
    """

    t_params = [IFG_LKSX, IFG_LKSY, IFG_CROP_OPT]
    xlooks, ylooks, crop = [params[k] for k in t_params]
    return xlooks, ylooks, crop


def original_ifg_paths(ifglist_path, obs_dir):
    """
    Returns sequence of paths to files in given ifglist file.

    Args:
        ifglist_path: Absolute path to interferogram file list.
        obs_dir: Absolute path to observations directory.

    Returns:
        list: List of full paths to interferogram files.
    """
    ifglist = parse_namelist(ifglist_path)
    return [os.path.join(obs_dir, p) for p in ifglist]


def coherence_paths_for(path: str, params: dict, tif=False) -> str:
    """
    Returns path to coherence file for given interferogram. Pattern matches
    based on epoch in filename.

    Example:
        '20151025-20160501_eqa_filt.cc'
        Date pair is the epoch.

    Args:
        path: Path to intergerogram to find coherence file for.
        params: Parameter dictionary.
        tif: Find converted tif if True (_cc.tif), else find .cc file.

    Returns:
        Path to coherence file.
    """
    _, filename = split(path)
    epoch = re.match(sixteen_digits_pattern, filename).group(0)
    if tif:
        coh_file_paths = [f.converted_path for f in params[COHERENCE_FILE_PATHS] if epoch in f.converted_path]
    else:
        coh_file_paths = [f.unwrapped_path for f in params[COHERENCE_FILE_PATHS] if epoch in f.unwrapped_path]

    if len(coh_file_paths) > 2:
        raise ConfigException(f"'{COH_FILE_DIR}': found more than one coherence "
                      f"file for '{path}'. There must be only one "
                      f"coherence file per interferogram. Found {coh_file_paths}.")
    return coh_file_paths[0]


def mlooked_path(path, looks, crop_out):
    """
    Adds suffix to ifg path, for creating a new path for multilooked files.

    :param str path: original interferogram path
    :param int looks: number of range looks applied
    :param int crop_out: crop option applied

    :return: multilooked file name
    :rtype: str
    """
    base, ext = splitext(path)
    return "{base}_{looks}rlks_{crop_out}cr{ext}".format(base=base, looks=looks, crop_out=crop_out, ext=ext)


def get_dest_paths(base_paths, crop, params, looks):
    """
    Determines the full path names for the destination multilooked files

    :param list base_paths: original interferogram paths
    :param int crop: Crop option applied
    :param dict params: Parameters dictionary
    :param int looks: number of range looks applied

    :return: full path names for destination files
    :rtype: list
    """

    dest_mlooked_ifgs = [mlooked_path(os.path.basename(q).split('.')[0] + '_'
                                      + os.path.basename(q).split('.')[1] +
                                      '.tif', looks=looks, crop_out=crop)
                         for q in base_paths]

    return [os.path.join(params[OUT_DIR], p) for p in dest_mlooked_ifgs]


# ==== PARAMETER VALIDATION ==== #

_PARAM_VALIDATION = {
    IFG_FILE_LIST: (
        lambda a: a is not None and os.path.exists(a),
        f"'{IFG_FILE_LIST}': file must be provided and must exist."
    ),
    DEM_FILE: (
        lambda a: a is not None and os.path.exists(a),
        f"'{DEM_FILE}': file must be provided and must exist."
    ),
    DEM_HEADER_FILE: (
        lambda a: a is not None and os.path.exists(a),
        f"'{DEM_HEADER_FILE}': file must be provided and must exist."
    ),
    OUT_DIR: (
        lambda a: a is not None,
        f"'{OBS_DIR}': directory must be provided."
    ),
    APS_INCIDENCE_MAP: (
        lambda a: os.path.exists(a) if a is not None else True,
        f"'{APS_INCIDENCE_MAP}': file must exist."
    ),
    APS_ELEVATION_MAP: (
        lambda a: os.path.exists(a) if a is not None else True,
        f"'{APS_ELEVATION_MAP}': file must exists."
    ),
    IFG_CROP_OPT: (
        lambda a: a in (1, 2, 3, 4),
        f"'{IFG_CROP_OPT}': must select option 1, 2, 3, or 4."
    ),
    IFG_LKSX: (
        lambda a: a >= 1,
        f"'{IFG_LKSX}': must be >= 1."
    ),
    IFG_LKSY: (
        lambda a: a >= 1,
        f"'{IFG_LKSY}': must be >= 1."
    ),
    NO_DATA_VALUE: (
        lambda a: True,
        "Any float value valid."
    ),
    COH_MASK: (
        lambda a: a in (0, 1),
        f"'{COH_MASK}': must select option 0 or 1."
    ),
    REFX: (
        lambda a: True,
        "Any float value valid."
    ),
    REFY: (
        lambda a: True,
        "Any float value valid."
    ),
    ORBITAL_FIT: (
        lambda a: a in (0, 1),
        f"'{ORBITAL_FIT}': must select option 0 or 1."
    ),
    LR_NSIG: (
        lambda a: 1 <= a <= 10,
        f"'{LR_NSIG}': must be between 1 and 10 (inclusive)."
    ),
    LR_PTHRESH: (
        lambda a: a >= 1,
        f"'{LR_PTHRESH}': must be >= 1"
    ),
    LR_MAXSIG: (
        lambda a: 0 <= a <= 1000,
        f"'{LR_MAXSIG}': must be between 0 and 1000 (inclusive)."
    ),
    APSEST: (
        lambda a: a in (0, 1),
        f"'{APSEST}': must select option 0 or 1."
    ),
    TIME_SERIES_CAL: (
        lambda a: a in (0, 1),
        f"'{TIME_SERIES_CAL}': must select option 0 or 1."
    ),
    PARALLEL: (
        lambda a: a in (0, 1),
        f"'{PARALLEL}': must select option 0 or 1."
    ),
    PROCESSES: (
        lambda a: a >= 1,
        f"'{PROCESSES}': must be >= 1."
    ),
    PROCESSOR: (
        lambda a: a in (0, 1, 2),
        f"'{PROCESSOR}': must select option 0 or 1."
    ),
    NAN_CONVERSION: (
        lambda a: a in (0, 1),
        f"'{NAN_CONVERSION}': must select option 0 or 1."
    ),
    NO_DATA_AVERAGING_THRESHOLD: (
        lambda a: True,
        "Any float value valid."),
}
"""dict: basic validation functions for compulsory parameters."""

_CUSTOM_CROP_VALIDATION = {
    IFG_XFIRST: (
        lambda a: a is not None,
        f"'{IFG_XFIRST}': must be provided."
    ),
    IFG_XLAST: (
        lambda a: a is not None,
        f"'{IFG_XLAST}': must be provided."
    ),
    IFG_YFIRST: (
        lambda a: a is not None,
        f"'{IFG_YFIRST}': must be provided."
    ),
    IFG_YLAST: (
        lambda a: a is not None,
        f"'{IFG_YLAST}': must be provided.."
    ),
}
"""dict: basic validation functions for custom cropping parameters."""

_GAMMA_VALIDATION = {
    HDR_FILE_LIST: (
        lambda a: a is not None and os.path.exists(a),
        f"'{HDR_FILE_LIST}': file must be provided and must exist."
    ),
}
"""dict: basic validation functions for gamma parameters."""

_COHERENCE_VALIDATION = {
    COH_THRESH: (
        lambda a: 0.0 <= a <= 1.0,
        f"'{COH_THRESH}': must be between 0.0 and 1.0 (inclusive)."
    ),
    COH_FILE_LIST: (
        lambda a: a is not None and not os.path.exists(a),
        f"'{COH_FILE_LIST}': if file is provided it must exist."
    ),
}
"""dict: basic validation functions for coherence parameters."""

_ORBITAL_FIT_VALIDATION = {
    ORBITAL_FIT_METHOD: (
        lambda a: a in (1, 2),
        f"'{ORBITAL_FIT_METHOD}': must select option 1 or 2."
    ),
    ORBITAL_FIT_DEGREE: (
        lambda a: a in (1, 2, 3),
        f"'{ORBITAL_FIT_DEGREE}': must select option 1, 2 or 3."
    ),
    ORBITAL_FIT_LOOKS_X: (
        lambda a: a >= 1,
        f"'{ORBITAL_FIT_LOOKS_X}': must be >= 1."
    ),
    ORBITAL_FIT_LOOKS_Y: (
        lambda a: a >= 1,
        f"'{ORBITAL_FIT_LOOKS_Y}': must be >= 1."
    ),
}
"""dict: basic validation fucntions for orbital error correction parameters."""

_APSEST_VALIDATION = {
    TLPF_METHOD: (
        lambda a: a in (1, 2, 3),
        f"'{TLPF_METHOD}': must select option 1, 2 or 3."
    ),
    TLPF_CUTOFF: (
        lambda a: a >= YEARS_PER_DAY, # 1 day in years
        f"'{TLPF_CUTOFF}': must be >= {YEARS_PER_DAY}."
    ),
    TLPF_PTHR: (
        lambda a: a >= 1,
        f"'{TLPF_PTHR}': must be >= 1."
    ),
    SLPF_METHOD: (
        lambda a: a in (1, 2),
        f"'{SLPF_METHOD}': must select option 1 or 2."
    ),
    SLPF_CUTOFF: (
        lambda a: a >= 0.001,
        f"'{SLPF_CUTOFF}': must be >= 0.001."
    ),
    SLPF_ORDER: (
        lambda a: 1 <= a <= 3,
        f"'{SLPF_ORDER}': must be between 1 and 3 (inclusive)."
    ),
    SLPF_NANFILL: (
        lambda a: a in (0, 1),
        f"'{SLPF_NANFILL}': must select option 0 or 1."
    ),
}
"""dict: basic validation functions for atmospheric correction parameters."""

_TIME_SERIES_VALIDATION = {
    TIME_SERIES_PTHRESH: (
        lambda a: a >= 1,
        f"'{TIME_SERIES_PTHRESH}': must be >= 1."
    ),
    TIME_SERIES_SM_FACTOR: (
        lambda a: -5.0 <= a <= 0,
        f"'{TIME_SERIES_SM_FACTOR}': must be between -5.0 and 0."
    ),
    TIME_SERIES_SM_ORDER: (
        lambda a: a in (1, 2),
        f"'{TIME_SERIES_SM_ORDER}': must select option 1 or 2."
    ),
    TIME_SERIES_METHOD: (
        lambda a: a in (1, 2),
        f"'{TIME_SERIES_METHOD}': must select option 1 or 2."
    ),
}
"""dict: basic vaidation functions for time series parameters."""

_REFERENCE_PIXEL_VALIDATION = {
    REFNX: (
        lambda a: 1 <= a <= 50,
        f"'{REFNX}': must be between 1 and 50 (inclusive)."
    ),
    REFNY: (
        lambda a: 1 <= a <= 50,
        f"'{REFNY}': must be between 1 and 50 (inclusive)."
    ),
    REF_CHIP_SIZE: (
        lambda a: 1 <= a <= 101 and a % 2 == 1,
        f"'{REF_CHIP_SIZE}': must be between 1 and 101 (inclusive) and be odd."
    ),
    REF_MIN_FRAC: (
        lambda a: 0.0 <= a <= 1.0,
        f"'{REF_MIN_FRAC}': must be between 0.0 and 1.0 "
        "(inclusive)."
    ),
    REF_EST_METHOD: (
        lambda a: a in (1, 2),
        f"'{REF_EST_METHOD}': must select option 1 or 2."
    ),
}


class ConfigException(Exception):
    """
    Default exception class for configuration errors.
    """
