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
in configuration files. Examples of PyRate configuration files are
provided in the configs/ directory
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


def get_config_params(path: str, validate: bool=True, step: str=CONV2TIF) -> Dict:
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
    params = _parse_conf_file(txt, validate, step)
    params[TMPDIR] = os.path.join(os.path.abspath(params[OUT_DIR]), 'tmpdir')

    return params


def _parse_conf_file(content, validate: bool=True, step: str=CONV2TIF) -> Dict:
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

    return _parse_pars(parameters, validate, step)


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


def _parse_pars(pars, validate: bool=True, step: str=CONV2TIF) -> Dict:
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

    if validate:
        validate_parameters(pars, step)
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


def get_ifg_paths(config_file, step=CONV2TIF):
    """
    Read the configuration file, extract interferogram file list and determine
    input and output interferogram path names.

    :param str config_file: Configuration file path

    :return: base_unw_paths: List of unwrapped inteferograms
    :return: dest_paths: List of multi-looked and cropped geotifs
    :return: params: Dictionary corresponding to the config file
    :rtype: list
    :rtype: list
    :rtype: dict
    """
    params = get_config_params(config_file, step=step)
    ifg_file_list = params.get(IFG_FILE_LIST)

    xlks, _, crop = transform_params(params)

    # base_unw_paths need to be geotiffed by conv2tif and multilooked by prepifg
    base_unw_paths = original_ifg_paths(ifg_file_list, params[OBS_DIR])

    # dest_paths are tifs that have been coherence masked (if enabled),
    #  cropped and multilooked

    if "tif" in base_unw_paths[0].split(".")[1]:

        dest_paths = get_dest_paths(base_unw_paths, crop, params, xlks)
        for i, dest_path in enumerate(dest_paths):
            dest_paths[i] = dest_path.replace("_tif", "")
    else:
        dest_paths = get_dest_paths(base_unw_paths, crop, params, xlks)

    return base_unw_paths, dest_paths, params

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
"""dict: basic validation functions for reference pixel search parameters."""

def convert_geographic_coordinate_to_pixel_value(refpx, refpy, transform):
    """
    Converts a lat/long coordinate to a pixel coordinate given the 
    geotransform of the image.

    Args:
        refpx: The longitude of the coordinate.
        refpx: The latitude of the coordinate.
        transform: The geotransform array of the image.

    Returns:
        Tuple of refpx, refpy in pixel values.
    """
    # transform = ifg.dataset.GetGeoTransform()

    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    refpx = int((refpx - xOrigin) / pixelWidth)
    refpy = int((yOrigin - refpy) / pixelHeight)

    return int(refpx), int(refpy)


def validate_parameters(pars: Dict, step: str=CONV2TIF):
    """
    Main validation function. Calls validation subfunctions and gathers
    some required variables for performing validation.

    Args:
        pars: The parameters dictionary.
        step: The current step of the PyRate workflow. Determines what 
            parameters to validate and what resources are available in
            terms of geotiffs.

    Raises:
        ConfigException: If errors occur during parameter validation.
    """
    is_GAMMA = pars[PROCESSOR] == GAMMA
    ifl = pars[IFG_FILE_LIST]

    # TODO: Call bounds checking functions based on the current step
    # Call basic bounds checking functions for all parameters.
    validate_compulsory_parameters(pars)
    validate_optional_parameters(pars)
    
    # Validate that provided filenames contain correct epochs and that
    #  the files exist.
    if is_GAMMA:
        validate_epochs(ifl, SIXTEEN_DIGIT_EPOCH_PAIR)
        validate_epochs(pars[HDR_FILE_LIST], EIGHT_DIGIT_EPOCH)
        validate_gamma_headers(ifl, pars[HDR_FILE_LIST])
    else:
        validate_epochs(ifl, TWELVE_DIGIT_EPOCH_PAIR)

    if step == CONV2TIF:   
        # Check that unwrapped interferograms exist.
        validate_ifgs(ifl)

    elif step == PREPIFG: 
        # Check full res geotiffs exist before continuing.
        validate_tifs_exist(ifl, pars[OUT_DIR])

        # Check the minimum number of epochs.
        n_epochs = 0
        with open(ifl, "r") as f:
            list_of_epoches = []
            for line in f.readlines():

                PTN = re.compile(r'\d{8}')  # match 8 digits for the dates

                # if ROI_PAC
                if 0 == pars["processor"]:
                    PTN = re.compile(r'\d{6}')  # match 8 digits for the dates

                epochs = PTN.findall(line.strip())
                list_of_epoches.extend(epochs)

        list_of_epoches = set(list_of_epoches)
        n_epochs = len(list_of_epoches)
        validate_minimum_epochs(n_epochs, MINIMUM_NUMBER_EPOCHS)

        # validate
        with open(ifl, "r") as f:
            # validate params for each geotiff
            for line in f.readlines():
                line = line.strip('\n')
                if line.endswith('.tif'):
                    tif_file_path = Path(pars[OUT_DIR]).joinpath(Path(line.strip()).name)
                else:
                    p = Path(line)
                    base, ext = p.stem, p.suffix[1:]
                    tif_file_path = Path(pars[OUT_DIR]).joinpath(base + "_" + ext + ".tif")

                if not tif_file_path.exists():
                    raise Exception("GeoTiff: " + tif_file_path.as_posix() + " not found.")

                raster = os.path.join(tif_file_path)
                gtif = gdal.Open(raster)

                latitudes = []
                longitudes = []

                for line in gdal.Info(gtif).split('\n'):
                    if "Size is" in line:
                        x_size, y_size = line.split("Size is")[1].split(",")
                        x_size, y_size = int(x_size.strip()), int(y_size.strip())

                    for line_tag in ["Upper Left", "Lower Left", "Upper Right", "Lower Right"]:
                        if line_tag in line:
                            latitude, longitude = line.split(")")[0].split("(")[1].split(",")
                            latitudes.append(float(latitude.strip()))
                            longitudes.append(float(longitude.strip()))

                # validate multi-look parameters
                if pars["ifglksx"] < 0 or pars["ifglksx"] > x_size:
                    raise Exception("Value of ifglksx: "+str(pars["ifglksx"])+" out of bounds: [0,"+str(x_size)+"]")
                if pars["ifglksy"] < 0 or pars["ifglksy"] > x_size:
                    raise Exception("Value of ifglksy: "+str(pars["ifglksy"])+" out of bounds: [0,"+str(x_size)+"]")

                # validate crop parameters
                if pars["ifgxfirst"] < min(latitudes) or pars["ifgxfirst"] > max(latitudes):
                    raise Exception(
                        "ifgxfirst: " + str(pars["ifgxfirst"]) + " not with in range {" + str(
                            min(latitudes)) + "," + str(max(latitudes)) + "}"
                    )

                if pars["ifgxlast"] < min(latitudes) or pars["ifgxlast"] > max(latitudes):
                    raise Exception(
                        "ifgxlast: " + str(pars["ifgxlast"]) + " not with in range {" + str(min(latitudes)) + "," + str(
                            max(latitudes)) + "}"
                    )

                if pars["ifgyfirst"] < min(longitudes) or pars["ifgyfirst"] > max(longitudes):
                    raise Exception(
                        "ifgyfirst: " + str(pars["ifgyfirst"]) + " not with in range {" + str(
                            min(longitudes)) + "," + str(max(longitudes)) + "}"
                    )

                if pars["ifgylast"] < min(longitudes) or pars["ifgylast"] > max(longitudes):
                    raise Exception(
                        "ifgylast: " + str(pars["ifgylast"]) + " not with in range {" + str(
                            min(longitudes)) + "," + str(max(longitudes)) + "}"
                    )

                del gtif  # manually close raster

        # Check coherence masking if enabled
        if pars[COH_MASK]:
            validate_epochs(pars[COH_FILE_LIST], SIXTEEN_DIGIT_EPOCH_PAIR)
            validate_coherence_files(ifl, pars)

    elif step == PROCESS:

        tmp_ifg_list = os.path.join(pars[OUT_DIR], "tmpIfgList")
        with open(tmp_ifg_list, "w") as fileHandler:
            for p in Path(pars[OUT_DIR]).rglob("*rlks_*cr.tif"):
                if "dem" not in str(p):
                    fileHandler.write(p.name + "\n")
        # Check that cropped/subsampled tifs exist.
        validate_prepifg_tifs_exist(tmp_ifg_list, pars[OUT_DIR], pars)

        # Check the minimum number of epochs.
        n_epochs, max_span = _get_temporal_info(tmp_ifg_list, pars[OBS_DIR])
        validate_minimum_epochs(n_epochs, MINIMUM_NUMBER_EPOCHS)

        # Get spatial information from tifs.

        extents, n_cols, n_rows, transform = _get_prepifg_info(tmp_ifg_list, pars[OBS_DIR], pars)

        # test if refx/y already set to default value of -1
        if pars[REFX] != -1 and pars[REFY] != -1:

            # Convert refx/refy from lat/long to pixel and validate...
            if (pars[REFX] <= 180) and (pars[REFX] >= -180) and (pars[REFY] >= -90) and (pars[REFY] <= 90):

                pars[REFX], pars[REFY] = convert_geographic_coordinate_to_pixel_value(pars[REFX], pars[REFY], transform)
                _logger.debug("converted pars[REFX], pars[REFY] to: " + str(pars[REFX]) + " " + str(pars[REFY]))

                if pars[REFX] < 0 or pars[REFX] > n_cols:
                    _logger.info("converted pars[REFX] out of range")
                    pars[REFX] = -1
                if pars[REFY] < 0 or pars[REFY] > n_rows:
                    _logger.info("converted pars[REFY] out of range")
                    pars[REFY] = -1

            # otherwise we need to search for the pixel so validate the search parameters.
            else:
                _logger.info("given pars[REFX] or pars[REFY] out of range")
                pars[REFX] = -1
                pars[REFY] = -1
                # TODO fix tests when below validation is included
                # validate_reference_pixel_search_windows(n_cols, n_rows, pars)

        validate_multilook_parameters(n_cols, n_rows, ORBITAL_FIT_LOOKS_X, ORBITAL_FIT_LOOKS_Y, pars)
        validate_slpf_cutoff(extents, pars)
        validate_epoch_thresholds(n_epochs, pars)
        validate_epoch_cutoff(max_span, TLPF_CUTOFF, pars)
        validate_obs_thresholds(tmp_ifg_list, pars)

    elif step == MERGE:
        validate_prepifg_tifs_exist(ifl, pars[OBS_DIR], pars)


def _raise_errors(errors: List[str]) -> bool:
    """
    Convenience function for raising an exception with errors.
    """
    if errors:
        errors.insert(0, "invalid parameters")
        raise ConfigException('\n'.join(errors))
    return True

def validate_compulsory_parameters(pars: Dict) -> Optional[bool]:
    """
    Calls the validators for compulsory (always used) parameters.

    Args:
        pars: The parameters dictionary.

    Returns:
        True if validation is successful.

    Raises:
        ConfigException: If validation fails.
    """
    errors = []
    for k in pars.keys():
        validator = _PARAM_VALIDATION.get(k)
        if validator is None:
            # _logger.debug(f"No validator implemented for '{k}'.")
            continue
        if not validator[0](pars[k]):
            errors.append(validator[1])

    return _raise_errors(errors)

def validate_optional_parameters(pars: Dict):
    """
    Calls the validators for optional parameters.

    Args:
        pars: The parameters dictionary.

    Returns:
        True if validation successful.

    Raises:
        ConfigException: If validation fails.
    """
    def validate(on: bool, validators: Dict, pars: Dict) -> List[str]:
        """
        Convenience method for calling validators.

        Args:
            on: Determines whether to call the validators.
            validators: A dictionary of validator functions.
            pars: Parameters dictionary.

        Returns:
            A list of errors.
        """
        errors = []
        if on:
            for k, validator in validators.items():
                if not validator[0](pars[k]):
                    errors.append(validator[1])
        return errors

    errors = []

    errors.extend(validate(pars[COH_MASK], _COHERENCE_VALIDATION, pars))
    errors.extend(validate(pars[APSEST], _APSEST_VALIDATION, pars))
    errors.extend(validate(pars[TIME_SERIES_CAL], _TIME_SERIES_VALIDATION, pars))
    errors.extend(validate(pars[ORBITAL_FIT], _ORBITAL_FIT_VALIDATION, pars))
    errors.extend(validate(pars[PROCESSOR] == GAMMA, _GAMMA_VALIDATION, pars))
    errors.extend(validate(pars[IFG_CROP_OPT] == 3, _CUSTOM_CROP_VALIDATION, pars))
    # errors.extend(validate(pars[REFX] > 0 and pars[REFY] > 0, _REFERENCE_PIXEL_VALIDATION, pars))

    return _raise_errors(errors)

def validate_minimum_epochs(n_epochs: int, min_epochs: int) -> Optional[bool]:
    """
    Validates the minimum number of epochs required for PyRate to produce
    good results.

    Args:
        n_epochs: The number of unique epochs in the collection of interferograms
            provided as input.
        min_epochs: The minimum number of epochs PyRate requires.

    Returns:
        True if there are enough epochs to satisfy minimum epochs.

    Raises:
        ConfigException: If there are not enough epochs to satisfy min epochs.
    """
    errors = []
    if n_epochs < min_epochs:
        errors.append(f"'{IFG_FILE_LIST}': total number of epochs given is {n_epochs}."
                      f" {min_epochs} or more unique epochs are required by PyRate.")
    _raise_errors(errors)

def validate_epochs(file_list: str, pattern: str) -> Optional[bool]:
    """
    Validate that user provided file names contain the correct pattern of
    epochs.

    Args:
        file_list: Path to the list of files to validate.
        pattern: Regex string for finding the epoch(s).

    Returns:
        True if all names in file list contain the epoch pattern *once*.

    Raises:
        ConfigException: If not all names in the file list don't contain
            the epoch pattern or contain it more than once.
    """
    errors = []
    PTN = re.compile(pattern)
    filenames = parse_namelist(file_list)
    for fn in filenames:
        epochs = PTN.findall(fn)
        if not epochs:
            errors.append(f"'{file_list}': {fn} does not contain an epoch of "
                          f"format {pattern}.")
        if len(epochs) > 1:
            errors.append(f"'{file_list}': {fn} does contains more than "
                          f"one epoch of {pattern}. There must be only one "
                          f"epoch in the filename.")

    return _raise_errors(errors)

def validate_epoch_cutoff(max_span: float, cutoff: str, pars: Dict) -> Optional[bool]:
    """
    Validate cutoff parameters that rely on the data timespan.

    Args:
        max_span: The maximum temporal span of the provided data in years.
        cutoff: The key of the cutoff parameter.
        pars: The parameters dictionary.

    Returns:
        True if the cutoff is less than the maximum data timespan.

    Raises:
        ConfigException: If the cutoff is greater than the max data timespan.
    """
    errors = []
    if pars[cutoff] > max_span:
        errors.append("'{cutoff}': must be less than max time span of "
                      "data in years ({max_span}).")
    return _raise_errors(errors)

def validate_prepifg_tifs_exist(ifg_file_list: str, obs_dir: str, pars: Dict) -> Optional[bool]:
    """
    Validates that cropped and multilooked interferograms exist in geotiff
    format.

    Args:
        ifg_file_list: Path to file containing interfergram file names.
        obs_dir: Path to observations directory.
        pars: Parameters dictionary.

    Returns:
        True if all interferograms exist in geotiff format.

    Raises:
        ConfigException: If not all intergerograms exist in geotiff format.
    """

    errors = []
    base_paths = [os.path.join(obs_dir, ifg) for ifg in parse_namelist(ifg_file_list)]
    ifg_paths = get_dest_paths(base_paths, pars[IFG_CROP_OPT], pars, pars[IFG_LKSX])
    for i, ifg_path in enumerate(ifg_paths):
        ifg_paths[i] = ifg_path.replace("_tif", "")

    for path in ifg_paths:
        if not os.path.exists(path):
            fname = os.path.split(path)[1]
            errors.append(f"'{IFG_FILE_LIST}': interferogram '{fname}' is "
                          f"required as a cropped and subsampled geotiff but "
                          f"could not be found. Make sure 'prepifg' has been "
                          f"run and ensure the '{IFG_LKSX}' and '{IFG_CROP_OPT}' "
                          f"parameters have not been changed since 'prepifg' was run.")

    return _raise_errors(errors)


def validate_tifs_exist(ifg_file_list: str, obs_dir: str) -> Optional[bool]:
    """
    Validates that provided interferograms exist in geotiff format.

    Args:
        ifg_file_list: Path to file containing interfergram file names.
        obs_dir: Path to observations directory.

    Returns:
        True if all interferograms exist in geotiff format.

    Raises:
        ConfigException: If not all intergerograms exist in geotiff format.
    """
    from pyrate.core.shared import output_tiff_filename

    errors = []
    ifgs = parse_namelist(ifg_file_list)
    ifg_paths = [os.path.join(obs_dir, ifg) for ifg in ifgs]
    gtiff_paths = [output_tiff_filename(f, obs_dir) for f in ifg_paths]
    for gtp in gtiff_paths:
        if not os.path.exists(gtp):
            fname = os.path.split(gtp)[1]
            errors.append(f"'{IFG_FILE_LIST}': interferogram '{fname}' is "
                          "required in geotiff format but no geotiff file "
                          "could be found.")

    return _raise_errors(errors)


def validate_ifgs(ifg_file_list: str) -> Optional[bool]:
    """
    Validates that provided interferograms exist.

    Args:
        ifg_file_list: Path to file containing interferogram file names.
        obs_dir: Path to observations directory.

    Returns:
        True if all interferograms exist.

    Raises:
        ConfigException: If not all interferograms exist.
    """
    errors = []
    ifgs = parse_namelist(ifg_file_list)
    for path in ifgs:
        if not os.path.exists(path):
            fname = os.path.split(path)[1]
            errors.append(f"'{IFG_FILE_LIST}': interferogram '{fname}' does not exist.")

    return _raise_errors(errors)


def validate_obs_thresholds(ifg_file_list: str, pars: Dict) -> Optional[bool]:
    """
    Validates parameters that specify an observations threshold.

    Args:
        ifg_file_list: Path to the file containing interferogram file names.
        pars: Parameters dictionary.

    Returns:
        True if there are enough interferograms to satisfy all observation thresholds.

    Raises:
        ConfigException: If there not enough interferograms to satisfy all observation
            thresholds.
    """
    def validate(n, p, k):
        thresh = p[k]
        if thresh > n:
            return [f"'{k}': not enough interferograms have been specified "
                    f"({n}) to satisfy threshold ({thresh})."]
        return []

    errors = []
    n_ifgs = len(list(parse_namelist(ifg_file_list)))
    errors.extend(validate(n_ifgs, pars, TIME_SERIES_PTHRESH))
    if pars[APSEST]:
        errors.extend(validate(n_ifgs, pars, TLPF_PTHR))

    return _raise_errors(errors)

def validate_epoch_thresholds(n_epochs: int, pars: Dict) -> Optional[bool]:
    """
    Validates threshold paramters that rely on the number of epochs
    available.

    Args:
        n_epochs: The number of unique epochs in the collection of interferograms
            provided as input.
        pars: Parameters dictionary.

    Returns:
        True if there are enough epochs to satisfy all epoch thresholds.

    Raises:
        ConfigException: If there are not enough epochs to satisfy all epoch thresholds.
    """
    errors = []
    thresh = pars[LR_PTHRESH]
    if n_epochs < thresh:
        errors.append(f"'{LR_PTHRESH}': not enough epochs have been specified "
                      f"({n_epochs}) to satisfy threshold ({thresh}).")

    return _raise_errors(errors)

def validate_crop_parameters(min_extents: Tuple[float, float, float, float],
                             pars: Dict) -> Optional[bool]:
    """
    Validates IFG crop parameters by ensuring user provided crop coordinates
    are within the *minimum* (intersecting) bounds of the interferogram stack.
    
    Args:
        min_extents: The minimum extents of the interferogram stack. The
            crop coordinates must be within these extents.
        pars: Parameters dictionary.

    Returns:
        True if validation is successful.

    Raises:
        ConfigException: If validation fails.
    """
    errors = []
    min_xmin, min_ymin, min_xmax, min_ymax = min_extents
    x_dim_string = f"(xmin: {min_xmin}, xmax: {min_xmax})"
    y_dim_string = f"(ymin: {min_ymin}, ymax: {min_ymax})"

    # Check crop coordinates within scene.
    def _validate_crop_coord(var_name, dim_min, dim_max, dim_string):
        if not dim_min < pars[var_name] < dim_max:
            return [f"'{var_name}': crop coordinate ({pars[var_name]}) "
                    f"is outside bounds of scene {dim_string}."]

        return []
    
    errors.extend(_validate_crop_coord(IFG_XFIRST, min_xmin, min_xmax, x_dim_string))
    errors.extend(_validate_crop_coord(IFG_YFIRST, min_ymin, min_ymax, y_dim_string))
    errors.extend(_validate_crop_coord(IFG_XLAST, min_xmin, min_xmax, x_dim_string))
    errors.extend(_validate_crop_coord(IFG_YLAST, min_ymin, min_ymax, y_dim_string))
    
    return _raise_errors(errors)


def validate_slpf_cutoff(extents: Tuple[float, float, float, float], 
                               pars: Dict) -> Optional[bool]:
    """
    Validate SLPF_CUTOFF is within the bounds of the scene being processed
    (after prepifg cropping/subsampling has occurred).

    Args:
        extents : Tuple of (xmin, xmax, ymin, ymax) describing the extents
            of the scene, after cropping & multisampling, in degrees.
        pars: Parameters dictionary.

    Returns:
        True if validation is successful.

    Raises:
        ConfigException: If validation fails.
    """
    xmin, ymin, xmax, ymax = extents
    # Check SLPF_CUTOFF within scene *in kilometeres*.
    DEG_TO_KM = 111.32 # km per degree
    x_extent = abs(xmin - xmax)
    y_extent = abs(ymin - ymax)
    x_extent_km = x_extent *  DEG_TO_KM
    y_extent_km = y_extent *  DEG_TO_KM
    if pars[SLPF_CUTOFF] > max(x_extent_km, y_extent_km):
        errors = [f"'{SLPF_CUTOFF}': cutoff is out of bounds, must be "
                  "less than max scene bound (in km) "
                  f"({max(x_extent_km, y_extent_km)})."]
    else:
        errors = []

    return _raise_errors(errors)



def validate_multilook_parameters(cols: int, rows: int, 
                                 xlks_name: str, ylks_name: str, 
                                 pars: Dict) -> Optional[bool]:
    """
    Validate multilook parameters by ensuring resulting resolution will be
    at least 1 pixel and less than the current resolution.

    Args:
        cols: Number of pixel columns (X) in the scene.
        rows: Number of pixel rows (Y) in the scene.
        xlks_name: Key for looking up the X looks factor.
        ylks_name: Key for looking up the Y looks factor.
        pars: Parameters dictionary.

    Returns:
        True if validation is successful.

    Raises:
        ConfigException: If validation fails.
    """
    errors = []

    def _validate_mlook(var_name, dim_val, dim_str):
        if pars[var_name] > dim_val:
            return [f"'{var_name}': the multilook factor ({pars[var_name]}) "
                    f"must be less than the number of pixels on the {dim_str} "
                    f"axis ({dim_val})."]
        return []

    errors.extend(_validate_mlook(xlks_name, cols, 'x'))
    errors.extend(_validate_mlook(ylks_name, rows, 'y'))
    return _raise_errors(errors)

def validate_reference_pixel_search_windows(n_cols: int, n_rows: int,
                                            pars: Dict) -> Optional[bool]:
    """
    Validates that the reference pixel search windows provided by user
    fit within the scene being processed.

    Args:
        n_cols: Number of pixel columns (X) in the scene.
        n_rows: Number of pixel rows (Y) in the scene.

    Returns:
        True if the scene can accomodate the search windows (no overlap).

    Raises:
        ConfigException: If the scene cannot accomodate the search windows.
    """
    from math import floor

    errors = []
    refnx = pars[REFNX]
    refny = pars[REFNY]
    chip_size = pars[REF_CHIP_SIZE]

    x_windows = floor(n_cols/chip_size)
    if refnx > x_windows:
        errors.append(f"'{REFNX}' & '{REF_CHIP_SIZE}': search windows do not "
                      f"fit in scene on X axis (number of columns: {n_cols}). "
                      f"Reduce {REF_CHIP_SIZE} or {REFNX} so that {REFNX} "
                      f"is less than or equal to (columns / {REF_CHIP_SIZE}).")
    y_windows = floor(n_rows/chip_size)
    if refny > y_windows:
        errors.append(f"'{REFNY}' & '{REF_CHIP_SIZE}': search windows do not "
                      f"fit in scene on Y axis (number of rows: {n_rows}). "
                      f"Reduce {REF_CHIP_SIZE} or {REFNY} so that {REFNY} "
                      f"is less than or equal to (rows / {REF_CHIP_SIZE}).")

    return _raise_errors(errors)


def validate_gamma_headers(ifg_file_list: str, slc_file_list: str) -> Optional[bool]:
    """
    Validates that a pair of GAMMA headers exist for each provided
    GAMMA interferogram.

    Args:
        ifg_file_list: Path to the file listing filenames of interferograms.
        slc_file_list: Path to the file listing filenames of GAMMA headers.
        slc_dir: Path to directory containing GAMMA headers.

    Returns:
        True if there are exactly 2 headers for each interogram (one for
        each epoch).

    Raises:
        ConfigException: If there are 0 or more than 2 matching headers
            for an interferogram.
    """
    from pyrate.core.gamma import get_header_paths
    errors = []

    for ifg in parse_namelist(ifg_file_list):
        headers = get_header_paths(ifg, slc_file_list)
        if len(headers) < 2:
            errors.append(f"'{SLC_DIR}': Headers not found for interferogram '{ifg}'.")

    return _raise_errors(errors)


def validate_coherence_files(ifg_file_list: str, pars: Dict) -> Optional[bool]:
    """
    Validates that there is a matching coherence file for each provided
    interferogram.

    Args:
        ifg_file_list: Path to file listing interferogram names.
        pars: The parameters dictionary.

    Returns:
        True if there is exactly 1 coherence file for each interferogram.

    Raises:
        ConfigException: If there are 0 or more than 1 matching coherence
            files for an interferogram.
    """
    errors = []

    for ifg in parse_namelist(ifg_file_list):
        paths = coherence_paths_for(ifg, pars)
        if not paths:
            errors.append(f"'{COH_FILE_DIR}': no coherence files found for intergerogram '{ifg}'.")
    return _raise_errors(errors)


def _get_temporal_info(ifg_file_list: str, obs_dir: str) -> Tuple:
    """
    Retrieves the number of unique epochs and maximum time span of the 
    dataset.

    Args:
        ifg_file_list: Path to file containing list of interferogram file names.
        obs_dir: Path to observations directory.
    
    Returns:
        Tuple containing the number of unique epochs and the maximum timespan.
    """
    from pyrate.core.algorithm import get_epochs
    from pyrate.core.shared import Ifg, output_tiff_filename

    ifg_paths = [os.path.join(obs_dir, ifg) for ifg in parse_namelist(ifg_file_list)]
    rasters = [Ifg(output_tiff_filename(f, obs_dir)) for f in ifg_paths]

    for r in rasters:
        if not r.is_open:
            r.open(readonly=True)

    # extents = xmin, ymin, xmax, ymax
    epoch_list = get_epochs(rasters)[0]
    n_epochs = len(epoch_list.dates)
    max_span = max(epoch_list.spans)

    for r in rasters:
        if not r.is_open:
            r.close()

    return n_epochs, max_span


def _get_prepifg_info(ifg_file_list: str, obs_dir: str, pars: Dict) -> Tuple:
    """
    Retrives spatial information from prepifg interferograms (images that
    have been cropped and multilooked).

    Args:
        ifg_file_list: Path to file containing list of interferogram file names.
        obs_dir: Path to observations directory.

    Returns:
    """
    from pyrate.core.shared import Ifg
    base_paths = [os.path.join(obs_dir, ifg) for ifg in parse_namelist(ifg_file_list)]
    ifg_paths = get_dest_paths(base_paths, pars[IFG_CROP_OPT], pars, pars[IFG_LKSX])

    for i, ifg_path in enumerate(ifg_paths):
        ifg_paths[i] = ifg_path.replace("_tif","")

    # This function assumes the stack of interferograms now have the same
    # bounds and resolution after going through prepifg.
    raster = Ifg(ifg_paths[0])
    if not raster.is_open:
        raster.open(readonly=True)

    n_cols = raster.ncols
    n_rows = raster.nrows
    extents = raster.x_first, raster.y_first, raster.x_last, raster.y_last
    transform = raster.dataset.GetGeoTransform()

    raster.close()

    return extents, n_cols, n_rows, transform


def _get_fullres_info(ifg_file_list: str, out_dir: str, crop_opts: Tuple) -> Tuple:
    """
    Retrieves spatial information from the provided interferograms.
    Requires the interferograms to exist in geotiff format.

    Args:
        ifg_file_list: Path to file containing list of interferogram file names.
        out_dir: Path to observations directory.

    Returns:
        Tuple containing extents (xmin, ymin, xmax, ymax), number of pixel
        columns, number of pixel rows, number of unique epochs and maximum
        time span of the data.
    """
    from pyrate.core.prepifg_helper import _min_bounds, _get_extents
    from pyrate.core.shared import Ifg, output_tiff_filename

    ifg_paths = parse_namelist(ifg_file_list)
    rasters = [Ifg(output_tiff_filename(f, out_dir)) for f in ifg_paths]

    for r in rasters:
        if not r.is_open:
            r.open(readonly=True)

    # extents = xmin, ymin, xmax, ymax
    min_extents = _min_bounds(rasters)
    post_crop_extents = _get_extents(rasters, crop_opts[0], user_exts=crop_opts[1])
    x_step = rasters[0].x_step
    y_step = rasters[0].y_step
    # Get the pixel bounds. Ifg/Raster objects do have 'ncols'/'nrows'
    # properties, but we'll calculate it off the extents we got above
    # because these take into account the chosen cropping option (until
    # the stack of interferograms is cropped it's not known what the
    # pixel dimensions will be).
    n_cols = abs(int(abs(post_crop_extents[0] - post_crop_extents[2]) / x_step))
    n_rows = abs(int(abs(post_crop_extents[1] - post_crop_extents[3]) / y_step))

    for r in rasters:
        r.close()

    return min_extents, n_cols, n_rows

def _crop_opts(params: Dict) -> Tuple:
    """
    Convenience function for getting crop options from parameters.
    """
    from pyrate.core.prepifg_helper import CustomExts

    crop_opt = params[IFG_CROP_OPT]
    if crop_opt == 3:
        xfirst = params[IFG_XFIRST]
        yfirst = params[IFG_YFIRST]
        xlast = params[IFG_XLAST]
        ylast = params[IFG_YLAST]
        return crop_opt, CustomExts(xfirst, yfirst, xlast, ylast)

    return crop_opt, None

class ConfigException(Exception):
    """
    Default exception class for configuration errors.
    """
