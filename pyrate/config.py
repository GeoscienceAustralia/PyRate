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
from __future__ import print_function
# TODO: add regex column to check if some values are within bounds? Potential
# problem with the checking being done in the middle of the runs, as bad values
# could cause crashes & destroying some of the results.
import os
from os.path import splitext
import warnings
from pyrate import compat
from pyrate import mpiops

# general constants
NO_MULTILOOKING = 1

# interp is automatically assigned; not needed in conf file
#TIME_SERIES_INTERP = 'tsinterp'

# constants for lookups
#: STR; Name of input interferogram list file
IFG_FILE_LIST = 'ifgfilelist'
#: STR; The name of the interferogram processor used (0==ROIPAC, 1==GAMMA)
PROCESSOR = 'processor'
#: STR; Name of directory containing input interferograms.
#: In the case of Python PyRate, these are the tif files,
#: Not the outputs from gamma or roipac.
OBS_DIR = 'obsdir'
#: STR; Name of directory for saving output products
OUT_DIR = 'outdir'
#: INT; Number of simulated datasets NOT CURRENTLY USED
#NUM_SETS = 'nsets'
#: STR; Directory containing simulated datasets NOT CURRENTLY USED
SIM_DIR = 'simdir'
#: STR; Name of Digital Elevation Model file used in
# constructing the interferograms
DEM_FILE = 'demfile'
#: STR; Name of the header for the DEM
DEM_HEADER_FILE = 'demHeaderFile'

# location of gamma slc files
SLC_DIR = 'slcFileDir'

#: STR; The projection of the input interferograms.
# When *PROCESSOR* == 0, either
#: this or *ROIPAC_RESOURCE_HEADER* must be provided.
INPUT_IFG_PROJECTION = 'projection'
#: STR; The resource header used for conferting ROIPAC interferograms.
# When *PROCESSOR* == 0, either
#: this or *INPUT_IFG_PROJECTION* must be provided.
ROIPAC_RESOURCE_HEADER = 'resourceHeader'
#: FLOAT; The no data value in the interferogram files.
NO_DATA_VALUE = 'noDataValue'
#: FLOAT; No data averaging threshold for prepifg
NO_DATA_AVERAGING_THRESHOLD = 'noDataAveragingThreshold'
#: BOOL (1/0); Use amplitude images NOT CURRENTLY USED
AMPLITUDE_FLAG = 'ampflag'
#: BOOL (1/0); Use baseline information NOT CURRENTLY USED
PERP_BASELINE_FLAG = 'basepflag'
#: BOOL (1/2/3); Re-project data from Line of sight, 1 = vertical,
# 2 = horizontal, 3 = no conversion
REPROJECTION = 'prjflag'
#: BOOL (0/1); Select MST algorithm, 0 = Matlab-Pirate algorithm, 1 = NetworkX
NETWORKX_OR_MATLAB_FLAG = 'networkx_or_matlab'
#: TODO; what does this parameter do?
NAN_CONVERSION = 'nan_conversion'

#: BOOL (1/2/3/4); Method for cropping interferograms,
# 1 = minimum overlapping area (intersection), 2 = maximum area (union),
# 3 = customised area, 4 = all ifgs already same size
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
#: FLOAT; Coordinate in x of reference pixel OR -1 = perform
# reference pixel search
REFX = 'refx'
#: FLOAT; Coordinate in y of reference pixel OR -1 = perform
# reference pixel search
REFY = 'refy'
#: INT; Number of reference pixel grid search nodes in x dimension
REFNX = "refnx"
#: INT; Number of reference pixel grid search nodes in y dimension
REFNY = "refny"
#: INT; Dimension of reference pixel search window
REF_CHIP_SIZE = 'refchipsize'
#: REAL; Minimum fraction of observations required in
# reference pixel search window for pixel to be a viable reference pixel
REF_MIN_FRAC = 'refminfrac'
# REFERENCE estimation method
REF_EST_METHOD = 'refest'

#atmospheric error correction parameter
APS_CORRECTION = 'apscorrect'
APS_METHOD = 'apsmethod'
APS_INCIDENCE_MAP = 'incidencemap'
APS_INCIDENCE_EXT = 'APS_INCIDENCE_EXT'
APS_ELEVATION_MAP = 'elevationmap'
APS_ELEVATION_EXT = 'APS_ELEVATION_EXT'

# orbital error correction/parameters
#: BOOL (1/0); Boolean flag controlling whether to apply orbital
# error correction
ORBITAL_FIT = 'orbfit'
#: BOOL (1/2); Method for orbital error correction, 1: ifg by ifg/independent,
#  2: epoch by epoch/network
ORBITAL_FIT_METHOD = 'orbfitmethod'
#: BOOL (1/2/3) Order of orbital error model, 1 = planar
# in x and y (2 parameter model, 2 = quadratic in x and y (5 parameter model),
# 3 = quadratic in x and cubic in y (part-cubic 6 parameter model)
ORBITAL_FIT_DEGREE = 'orbfitdegrees'
#: INT; Multi look factor for orbital error calculation in x dimension
ORBITAL_FIT_LOOKS_X = 'orbfitlksx'
#: INT; Multi look factor for orbital error calculation in y dimension
ORBITAL_FIT_LOOKS_Y = 'orbfitlksy'
# ORBITAL_FIT_orbrefest:	 1 BOOLEAN (1/0) # remove reference phase
# ORBITAL_FIT_ orbmaskflag:   1 BOOLEAN (1/0) # mask some patches
# for orbital correction

# Linear rate/stacking parameters
#: REAL; Threshold ratio between 'model minus observation'
#  residuals and a-priori observation standard deviations for
# linear rate estimate acceptance (otherwise remove furthest
# outlier and re-iterate)
LR_NSIG = 'nsig'
#: INT; Number of required input observations per pixel for the
# linear rate inversion
LR_PTHRESH = 'pthr'
#: REAL; Maximum allowable standard error for pixels in linear rate inversion.
LR_MAXSIG = 'maxsig'

#: BOOL (1/0); Do Time series calculation
TIME_SERIES_CAL = 'tscal'
#: INT (1/2); Method for time series inversion (1: Laplacian Smoothing; 2: SVD)
TIME_SERIES_METHOD = 'tsmethod'
#: INT; Number of required input observations
# per pixel for time series inversion
TIME_SERIES_PTHRESH = 'ts_pthr'
#: INT (1/2); Order of Laplacian smoothing operator, first or
# second order NOT CURRENTLY USED
TIME_SERIES_SM_ORDER = 'smorder'
#: REAL; Laplacian smoothing factor (0: calculate & plot L-curve;
# others: using the specific smoothing factor 10**smfactor) NOT CURRENTLY USED
TIME_SERIES_SM_FACTOR = 'smfactor'

# MULTIPROCESSING parameters
PARALLEL = 'parallel'
PROCESSES = 'processes'

# Luigi parameter
LUIGI = 'use_luigi'

# pickle output or not
PICKLE = 'pickle'

# ORBITAL ERROR correction constants
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2

PLANAR = 'PLANAR'
QUADRATIC = 'QUADRATIC'
PART_CUBIC = 'PART_CUBIC'


def degree_conv(deg):
    """
    Convenience: convert numerical degree to human readable string
    """

    degree = int(deg)
    if degree == 1:
        return PLANAR
    if degree == 2:
        return QUADRATIC
    if degree == 3:
        return PART_CUBIC
    raise ValueError("Orbital fit polynomial degree option not recognised")


def method_conv(meth):
    """
    Convenience: convert numerical method to human readable string
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
    #PERP_BASELINE_FLAG : (bool, True),
    #AMPLITUDE_FLAG : (bool, False),
    #NUM_SETS : (int, 1),
    REPROJECTION : (int, 3),
    IFG_CROP_OPT : (int, None), # TODO: default to ALREADY_SAME_SIZE?
    IFG_LKSX : (int, NO_MULTILOOKING),
    IFG_LKSY : (int, NO_MULTILOOKING),
    IFG_XFIRST : (float, None),
    IFG_XLAST : (float, None),
    IFG_YFIRST : (float, None),
    IFG_YLAST : (float, None),
    NO_DATA_VALUE: (float, 0.0),

    REFX: (int, -1),
    REFY: (int, -1),
    REFNX: (int, 5),  # was 50 in original Pirate code
    REFNY: (int, 5),  # was 50 in original Pirate code
    REF_CHIP_SIZE: (int, 3),  # defaults to 21 in orig
    REF_MIN_FRAC: (float, 0.8),  # uses Pirate default
    REF_EST_METHOD: (int, 1),  # ref phase estimation method

    #ORBITAL_FIT : (bool, False),
    ORBITAL_FIT: (int, 0),
    ORBITAL_FIT_METHOD: (method_conv, NETWORK_METHOD),
    ORBITAL_FIT_DEGREE: (degree_conv, QUADRATIC),
    ORBITAL_FIT_LOOKS_X: (int, NO_MULTILOOKING),
    ORBITAL_FIT_LOOKS_Y: (int, NO_MULTILOOKING),

    LR_NSIG : (int, 3), # Pirate default
    LR_PTHRESH : (int, 20), # should be based on nepochs since not every project may have 20 epochs
    LR_MAXSIG : (int, 2), # Pirate default

    #TIME_SERIES_CAL : (bool, False),
    TIME_SERIES_CAL: (int, 0),
    TIME_SERIES_PTHRESH: (int, None),
    TIME_SERIES_SM_FACTOR: (float, None),
    TIME_SERIES_SM_ORDER: (int, None),
    TIME_SERIES_METHOD: (int, 2), # Default to SVD method

    PARALLEL: (int, None),
    PROCESSES: (int, 8),

    PROCESSOR: (int, None),
    NETWORKX_OR_MATLAB_FLAG: (int, 1), # Default to NetworkX

    LUIGI: (int, 0),
    NAN_CONVERSION: (int, 0),
    NO_DATA_AVERAGING_THRESHOLD: (float, 0.0),
    APS_CORRECTION: (int, 0),
    APS_METHOD: (int, 1)
    }
    #TIME_SERIES_INTERP : (bool, False) # Automatically assigned in code


PATHS = [OBS_DIR, IFG_FILE_LIST, DEM_FILE,
         DEM_HEADER_FILE, OUT_DIR,
         ROIPAC_RESOURCE_HEADER,
         SLC_DIR,
         APS_INCIDENCE_MAP,
         APS_ELEVATION_MAP]

INT_KEYS = [APS_CORRECTION, APS_METHOD]


def get_config_params(path):
    """
    Returns a dict for the key:value pairs from the .conf file
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
    if mpiops.size > 1 and params[LUIGI] == 1:
        raise ConfigException('LUIGI with MPI not supported. Please '
                              'turn off LUIGI in config file or '
                              'use LUIGI without MPI')
    return params


def _parse_conf_file(content):
    """
    Parser for converting text content into a dict of parameters
    """

    def is_valid(line):
        """
        check if line is  not empty or has % or #
        """
        return line != "" and line[0] not in "%#"

    lines = [ln.split() for ln in content.split('\n') if is_valid(ln)]

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

    parameters = handle_pyaps_parameters(parameters)

    if not parameters:
        raise ConfigException('Cannot parse any parameters from config file')

    return _parse_pars(parameters)


def handle_pyaps_parameters(params):
    """
    Parameters
    ----------
    params: dict

    Returns
    ------
    params: dict
    """
    params[APS_INCIDENCE_EXT] = None
    params[APS_ELEVATION_EXT] = None

    if compat.PyAPS_INSTALLED:
        # define APS_INCIDENCE_EXT for gamma prepifg
        if ((params[APS_INCIDENCE_MAP] is not None) and
                (params[APS_ELEVATION_MAP] is not None)):
            warnings.warn('Both incidence and elevation map supplied. '
                          'Using the incidence map and ignoring elevation map')

        if (int(params[APS_CORRECTION]) and
                (int(params[APS_METHOD]) == 2) and
                ((params[APS_INCIDENCE_MAP] is None) and
                 (params[APS_ELEVATION_MAP] is None))):
            raise ConfigException('When doing APS correction using method 2,'
                                  'the incidence/elevation map method,'
                                  'one of incidence or elevation map must be '
                                  'provided')

    if params[APS_INCIDENCE_MAP] is not None:
        params[APS_INCIDENCE_EXT] = \
            os.path.basename(params[APS_INCIDENCE_MAP]).split('.')[-1]
        params[APS_ELEVATION_MAP] = None
        params[APS_ELEVATION_EXT] = None
        return params

    # define APS_ELEVATON_EXT for gamma prepifg
    if params[APS_ELEVATION_MAP] is not None:
        params[APS_ELEVATION_EXT] = \
            os.path.basename(params[APS_ELEVATION_MAP]).split('.')[-1]

    return params


def _parse_pars(pars):
    """
    Parses and converts config file params from text
    """
    for k in PARAM_CONVERSION:
        if k in pars:
            conversion_func = PARAM_CONVERSION[k][0]
            if conversion_func:
                pars[k] = conversion_func(pars[k])
        else:
            # revert empty options to default value
            if k in PARAM_CONVERSION:
                pars[k] = PARAM_CONVERSION[k][1]
    return pars


def parse_namelist(nml):
    """
    Parses name list file into array of paths
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
    takes a param object and write the config file. Reverse of get_conf_params
    :param params: params dictionary

    """
    with open(output_conf_file, 'w') as f:
        for k, v in params.items():
            if k == ORBITAL_FIT_DEGREE:
                v = reverse_degree_conv(v)
            if v is not None:
                f.write(''.join([k, ':\t', str(v), '\n']))
            else:
                f.write(''.join([k, ':\t', '', '\n']))


def reverse_degree_conv(v):
    """
    Convenience: convert numerical degree to human readable string
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


def transform_params(params):
    """
    Returns subset of config parameters for cropping and multilooking.
    """

    t_params = [IFG_LKSX, IFG_LKSY, IFG_CROP_OPT]
    xlooks, ylooks, crop = [params[k] for k in t_params]
    return xlooks, ylooks, crop


def original_ifg_paths(ifglist_path):
    """
    Returns sequence of paths to files in given ifglist file.
    """

    basedir = os.path.dirname(ifglist_path)
    ifglist = parse_namelist(ifglist_path)
    return [os.path.join(basedir, p) for p in ifglist]


def mlooked_path(path, looks, crop_out):
    """
    Adds suffix to path, for creating a new path for mlooked files.
    """
    base, ext = splitext(path)
    return "{base}_{looks}rlks_{crop_out}cr{ext}".format(
        base=base, looks=looks, crop_out=crop_out, ext=ext)


def get_dest_paths(base_paths, crop, params, looks):
    """
    Parameters
    ----------
    base_paths: list
        unwrapped interferrogram paths
    crop: int
        crop method to use
    params: dict
        params dict corresponding to config file used
    looks: int
        multilooking

    Return
    ------
        list of output geotifs
    """
    dest_mlooked_ifgs = [mlooked_path(os.path.basename(q).split('.')[0] + '_'
                                      + os.path.basename(q).split('.')[1] +
                                      '.tif', looks=looks, crop_out=crop)
                         for q in base_paths]

    return [os.path.join(params[OUT_DIR], p) for p in dest_mlooked_ifgs]


def get_ifg_paths(config_file):
    """
    Parameters
    ----------
    config_file: str
        config file path

    Return
    ------
    base_unw_paths: list
        list of unwrapped inteferrograms
    dest_paths: list
        list of multilooked and cropped geotifs
    pars: dict
        dictionary corresponding to the config file.
    """
    pars = get_config_params(config_file)
    ifg_file_list = pars.get(IFG_FILE_LIST)
    pars[IFG_FILE_LIST] = ifg_file_list

    if ifg_file_list is None:
        emsg = 'Error {code}: Interferogram list file name not provided ' \
               'or does not exist'.format(code=2)
        raise IOError(2, emsg)
    xlks, _, crop = transform_params(pars)

    # base_unw_paths need to be geotiffed and multilooked by run_prepifg
    base_unw_paths = original_ifg_paths(ifg_file_list)

    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = get_dest_paths(base_unw_paths, crop, pars, xlks)

    return base_unw_paths, dest_paths, pars
