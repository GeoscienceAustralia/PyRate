#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
# pylint: disable= invalid-name
import os
from os.path import splitext
import warnings
import pickle
from pyrate import mpiops

# TODO: add regex column to check if some values are within bounds? Potential
# problem with the checking being done in the middle of the runs, as bad values
# could cause crashes & destroying some of the results.

# general constants
NO_MULTILOOKING = 1

# constants for lookups
#: STR; Name of input interferogram list file
IFG_FILE_LIST = 'ifgfilelist'
#: BOOL (0/1); The interferogram processor used (0==ROIPAC, 1==GAMMA)
PROCESSOR = 'processor'
#: STR; Name of directory containing input interferograms.
OBS_DIR = 'obsdir'
#: STR; Name of directory for saving output products
OUT_DIR = 'outdir'
#: STR; Name of Digital Elevation Model file
DEM_FILE = 'demfile'
#: STR; Name of the header for the DEM
DEM_HEADER_FILE = 'demHeaderFile'
#: STR; Name of directory containing GAMMA SLC parameter files
SLC_DIR = 'slcFileDir'

#: STR; The projection of the input interferograms.
INPUT_IFG_PROJECTION = 'projection'
#: FLOAT; The no data value in the interferogram files.
NO_DATA_VALUE = 'noDataValue'
#: FLOAT; No data averaging threshold for prepifg
NO_DATA_AVERAGING_THRESHOLD = 'noDataAveragingThreshold'
#: BOOL (1/2/3); Re-project data from Line of sight, 1 = vertical,
# 2 = horizontal, 3 = no conversion
REPROJECTION = 'prjflag' # NOT CURRENTLY USED
#: BOOL (0/1); Select MST algorithm, 0 = Matlab Pirate algorithm, 1 = NetworkX
NETWORKX_OR_MATLAB_FLAG = 'networkx_or_matlab'
#: BOOL (0/1): Convert no data values to Nan
NAN_CONVERSION = 'nan_conversion'

# Prepifg parameters
#: BOOL (1/2/3/4); Method for cropping interferograms, 1 = minimum overlapping area (intersection), 2 = maximum area (union), 3 = customised area, 4 = all ifgs already same size
IFG_CROP_OPT = 'ifgcropopt'
#: INT; Multi look factor for interferogram preparation in x dimension
IFG_LKSX = 'ifglksx'
#: INT; Multi look factor for interferogram preparation in y dimension
IFG_LKSY = 'ifglksy'
#: REAL; Minimum longitude for cropping with method 3
IFG_XFIRST = 'ifgxfirst'
#: REAL; Maximum longitude for cropping with method 3
IFG_XLAST = 'ifgxlast'
#: REAL; Minimum latitude for cropping with method 3
IFG_YFIRST = 'ifgyfirst'
#: REAL; Maximum latitude for cropping with method 3
IFG_YLAST = 'ifgylast'

# reference pixel parameters
#: INT; Coordinate in x of reference pixel OR -1 = perform search
REFX = 'refx'
#: INT; Coordinate in y of reference pixel OR -1 = perform search
REFY = 'refy'
#: INT; Number of reference pixel grid search nodes in x dimension
REFNX = "refnx"
#: INT; Number of reference pixel grid search nodes in y dimension
REFNY = "refny"
#: INT; Dimension of reference pixel search window
REF_CHIP_SIZE = 'refchipsize'
#: REAL; Minimum fraction of observations required in search window for pixel to be a viable reference pixel
REF_MIN_FRAC = 'refminfrac'
#: BOOL (1/2); Reference phase estimation method
REF_EST_METHOD = 'refest'

#atmospheric error correction parameters NOT CURRENTLY USED
APS_CORRECTION = 'apscorrect'
APS_METHOD = 'apsmethod'
APS_INCIDENCE_MAP = 'incidencemap'
APS_INCIDENCE_EXT = 'APS_INCIDENCE_EXT'
APS_ELEVATION_MAP = 'elevationmap'
APS_ELEVATION_EXT = 'APS_ELEVATION_EXT'

# orbital error correction/parameters
#: BOOL (1/0); Flag controlling whether to apply orbital error correction
ORBITAL_FIT = 'orbfit'
#: BOOL (1/2); Method for orbital error correction, 1: independent, 2: network
ORBITAL_FIT_METHOD = 'orbfitmethod'
#: BOOL (1/2/3) Order of orbital error model, 1 = planar in x and y (2 parameter model, 2 = quadratic in x and y (5 parameter model), 3 = quadratic in x and cubic in y (part-cubic 6 parameter model)
ORBITAL_FIT_DEGREE = 'orbfitdegrees'
#: INT; Multi look factor for orbital error calculation in x dimension
ORBITAL_FIT_LOOKS_X = 'orbfitlksx'
#: INT; Multi look factor for orbital error calculation in y dimension
ORBITAL_FIT_LOOKS_Y = 'orbfitlksy'

# Linear rate/stacking parameters
#: REAL; Threshold ratio between 'model minus observation' residuals and a-priori observation standard deviations for linear rate estimate acceptance (otherwise remove furthest outlier and re-iterate)
LR_NSIG = 'nsig'
#: INT; Number of required observations per pixel for the linear rate inversion
LR_PTHRESH = 'pthr'
#: REAL; Maximum allowable standard error for pixels in linear rate inversion.
LR_MAXSIG = 'maxsig'

# atmospheric delay errors fitting parameters NOT CURRENTLY USED
# atmfitmethod = 1: interferogram by interferogram; atmfitmethod = 2, epoch by epoch
#ATM_FIT = 'atmfit'
#ATM_FIT_METHOD = 'atmfitmethod'

#: BOOL (0/1) Do spatio-temporal filter
APSEST = 'apsest'

# temporal low-pass filter parameters
TLPF_METHOD = 'tlpfmethod'
TLPF_CUTOFF = 'tlpfcutoff'
TLPF_PTHR = 'tlpfpthr'

# spatially correlated noise low-pass filter parameters
SLPF_METHOD = 'slpfmethod'
SLPF_CUTOFF = 'slpfcutoff'
SLPF_ORDER = 'slpforder'
SLPF_NANFILL = 'slpnanfill'
SLPF_NANFILL_METHOD = 'slpnanfill_method'

# Time series parameters
#: BOOL (1/0); Do Time series calculation
TIME_SERIES_CAL = 'tscal'
#: INT (1/2); Method for time series inversion (1: Laplacian Smoothing; 2: SVD)
TIME_SERIES_METHOD = 'tsmethod'
#: INT; Number of required input observations per pixel for time series inversion
TIME_SERIES_PTHRESH = 'ts_pthr'
#: INT (1/2); Order of Laplacian smoothing operator, first or # second order
TIME_SERIES_SM_ORDER = 'smorder'
#: REAL; Laplacian smoothing factor (values used is 10**smfactor)
TIME_SERIES_SM_FACTOR = 'smfactor'
# tsinterp is automatically assigned in the code; not needed in conf file
#TIME_SERIES_INTERP = 'tsinterp'

#: BOOL (0/1/2); Use parallelisation/Multi-threading
PARALLEL = 'parallel'
#: INT; Number of processes for multi-threading
PROCESSES = 'processes'

# Orbital error correction constants for conversion to readable flags
INDEPENDENT_METHOD = 'INDEPENDENT'
NETWORK_METHOD = 'NETWORK'
PLANAR = 'PLANAR'
QUADRATIC = 'QUADRATIC'
PART_CUBIC = 'PART_CUBIC'

# dir for temp files
TMPDIR = 'tmpdir'


def _orb_degree_conv(deg):
    """
    Convenience: convert numerical degree flag to human readable string
    """
    degree = int(deg)
    if degree == 1:
        return PLANAR
    if degree == 2:
        return QUADRATIC
    if degree == 3:
        return PART_CUBIC
    raise ValueError("Orbital fit polynomial degree option not recognised")


def _orb_method_conv(meth):
    """
    Convenience: convert numerical method flag to human readable string
    """
    method = int(meth)
    if method == 1:
        return INDEPENDENT_METHOD
    if method == 2:
        return NETWORK_METHOD
    raise ValueError("Orbital fit method not recognised")

# Lookup to help convert args to correct type/defaults
# format is	key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = {
    REPROJECTION : (int, 3), # Default no conversion, CONVERSION NOT IMPLEMENTED
    IFG_CROP_OPT : (int, 1), # default to area 'intersection' option
    IFG_LKSX : (int, NO_MULTILOOKING),
    IFG_LKSY : (int, NO_MULTILOOKING),
    IFG_XFIRST : (float, None),
    IFG_XLAST : (float, None),
    IFG_YFIRST : (float, None),
    IFG_YLAST : (float, None),
    NO_DATA_VALUE: (float, 0.0),

    REFX: (int, -1),
    REFY: (int, -1),
    REFNX: (int, 50),
    REFNY: (int, 50),
    REF_CHIP_SIZE: (int, 21),
    REF_MIN_FRAC: (float, 0.8),
    REF_EST_METHOD: (int, 1),  # default to average of whole image

    ORBITAL_FIT: (int, 0),
    ORBITAL_FIT_METHOD: (_orb_method_conv, NETWORK_METHOD),
    ORBITAL_FIT_DEGREE: (_orb_degree_conv, QUADRATIC),
    ORBITAL_FIT_LOOKS_X: (int, NO_MULTILOOKING),
    ORBITAL_FIT_LOOKS_Y: (int, NO_MULTILOOKING),

    LR_NSIG: (int, 3),
    # pixel thresh based on nepochs? not every project may have 20 epochs
    LR_PTHRESH: (int, 20),
    LR_MAXSIG: (int, 2),

    #ATM_FIT: (int, 0), NOT CURRENTLY USED
    #ATM_FIT_METHOD: (int, 2),

    APSEST: (int, 0),
    TLPF_METHOD: (int, 1),
    TLPF_CUTOFF: (float, 0.0),
    TLPF_PTHR: (int, 1),

    SLPF_METHOD: (int, 1),
    SLPF_CUTOFF: (float, 0.0),
    SLPF_ORDER: (int, 1),
    SLPF_NANFILL: (int, 0),

    TIME_SERIES_CAL: (int, 0),
    # pixel thresh based on nepochs? not every project may have 20 epochs
    TIME_SERIES_PTHRESH: (int, 20),
    TIME_SERIES_SM_FACTOR: (float, None),
    TIME_SERIES_SM_ORDER: (int, None),
    TIME_SERIES_METHOD: (int, 2),  # Default to SVD method

    PARALLEL: (int, 0),
    PROCESSES: (int, 8),
    PROCESSOR: (int, None),
    NETWORKX_OR_MATLAB_FLAG: (int, 1), # Default to NetworkX
    NAN_CONVERSION: (int, 0),
    NO_DATA_AVERAGING_THRESHOLD: (float, 0.0),
    APS_CORRECTION: (int, 0),
    APS_METHOD: (int, 1)
    }


PATHS = [OBS_DIR, IFG_FILE_LIST, DEM_FILE,
         DEM_HEADER_FILE, OUT_DIR,
         SLC_DIR,
         APS_INCIDENCE_MAP,
         APS_ELEVATION_MAP]

INT_KEYS = [APS_CORRECTION, APS_METHOD]


def get_config_params(path):
    """
    Returns a dictionary of key:value parameter pairs from the
    configuration file

    :param str path: path of config file

    :return: params
    :rtype: dict
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


def _parse_conf_file(content):
    """
    Parser for converting text content into a dictionary of parameters
    """

    def _is_valid(line):
        """
        Check if line is not empty or has % or #
        """
        return line != "" and line[0] not in "%#"

    lines = [ln.split() for ln in content.split('\n') if _is_valid(ln)]

    # convert "field:   value" lines to [field, value]
    kvpair = [(e[0].rstrip(":"), e[1]) for e in lines if len(e) == 2] \
        + [(e[0].rstrip(":"), None) for e in lines if len(e) == 1]
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
    Function to check if requirements for weather model correction are given
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
        params[APS_ELEVATION_EXT] = os.path.basename(params[APS_ELEVATION_MAP]).split('.')[-1]

    return params


def _parse_pars(pars):
    """
    Parses and converts config file params from text
    """
    for k in PARAM_CONVERSION:
        if k in pars:
            # if option value is blank/missing revert to default
            if pars[k] is None:
                pars[k] = PARAM_CONVERSION[k][1]
            conversion_func = PARAM_CONVERSION[k][0]
            if conversion_func:
                pars[k] = conversion_func(pars[k])
        else:
            # revert missing options to default value
            if k in PARAM_CONVERSION:
                pars[k] = PARAM_CONVERSION[k][1]
    return pars


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


class ConfigException(Exception):
    """
    Default exception class for configuration errors.
    """
    pass


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
            if k == ORBITAL_FIT_DEGREE:
                v = _reverse_orb_degree_conv(v)
            if k == ORBITAL_FIT_METHOD:
                v = _reverse_orb_method_conv(v)
            if v is not None:
                f.write(''.join([k, ':\t', str(v), '\n']))
            else:
                f.write(''.join([k, ':\t', '', '\n']))


def _reverse_orb_degree_conv(v):
    """
    Convenience: convert degree to integer for config file
    """
    if v == PLANAR:
        return 1
    if v == QUADRATIC:
        return 2
    if v == PART_CUBIC:
        return 3
    else:
        raise ValueError(
            "Orbital fit polynomial degree option not recognised")


def _reverse_orb_method_conv(v):
    """
    Convenience: convert method to integer for config file
    """
    if v == INDEPENDENT_METHOD:
        return 1
    if v == NETWORK_METHOD:
        return 2
    else:
        raise ValueError(
            "Orbital fit method option not recognised")


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


def original_ifg_paths(ifglist_path):
    """
    Returns sequence of paths to files in given ifglist file.

    :param str ifglist_path: full path to interferogram file list

    :return: full path to ifg files
    :rtype: list
    """

    basedir = os.path.dirname(ifglist_path)
    ifglist = parse_namelist(ifglist_path)
    return [os.path.join(basedir, p) for p in ifglist]


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
    return "{base}_{looks}rlks_{crop_out}cr{ext}".format(
        base=base, looks=looks, crop_out=crop_out, ext=ext)


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


def get_ifg_paths(config_file):
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
    params = get_config_params(config_file)
    ifg_file_list = params.get(IFG_FILE_LIST)
    params[IFG_FILE_LIST] = ifg_file_list

    if ifg_file_list is None:
        emsg = 'Error {code}: Interferogram list file name not provided ' \
               'or does not exist'.format(code=2)
        raise IOError(2, emsg)
    xlks, _, crop = transform_params(params)

    # base_unw_paths need to be geotiffed and multilooked by run_prepifg
    base_unw_paths = original_ifg_paths(ifg_file_list)

    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = get_dest_paths(base_unw_paths, crop, params, xlks)

    return base_unw_paths, dest_paths, params
