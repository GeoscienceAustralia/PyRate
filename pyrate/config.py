"""
Utilities to parse PyRate configuration files. Includes numerous general PyRate
constants relating to options in configuration files.

An example PyRate configuration file is as follows::

    obsdir:       obs/
    ifgfilelist:  obs/ifms.list
    demfile:      dem/sydney_trimmed.tif
    basepflag:    0
    outdir:       out/

    ifgcropopt:   4
    ifglksx:      1
    ifglksy:      1
    ifgxfirst:    152.0
    ifgxlast:     152.2
    ifgyfirst:    -4.2
    ifgylast:     -4.3 

    orbfit:        0
    orbfitmethod:  1
    orbfitdegrees: 1
    orbfitlksx:    5
    orbfitlksy:    5

    refx:          12
    refy:          38
    refnx:         5
    refny:         5
    refchipsize:   5
    refminfrac:    0.8

    tscal:         1
    ts_pthr:       10
    ts_interp:     1

    nsig:          3
    pthr:          20
    maxsig:        2

Created on 17/09/2012

.. codeauthor:: Ben Davies
"""

# TODO: add regex column to check if some values are within bounds? Potential
# problem with the checking being done in the middle of the runs, as bad values
# could cause crashes & destroying some of the results.

import pyrate.orbital as orbital


# general constants
NO_MULTILOOKING = 1

# constants for lookups
NUMBER_OF_SETS = 'nsets'
#: STR; Name of input interferogram file
IFG_FILE_LIST = 'ifgfilelist'
#: STR; Name of directory containing input interferograms
OBS_DIR = 'obsdir'
#: STR; Name of directory for saving output products
OUT_DIR = 'outdir'
#: INT; Number of simulated datasets NOT CURRENTLY USED
NUM_SETS = 'nsets'
#: STR; Directory containing simulated datasets NOT CURRENTLY USED
SIM_DIR = 'simdir'
#: STR; Name of Digital Elevation Model file
DEM_FILE = 'demfile'
#: BOOL (1/0); Use amplitude images NOT CURRENTLY USED 
AMPLITUDE_FLAG = 'ampflag'
#: BOOL (1/0); Use baseline information NOT CURRENTLY USED
PERP_BASELINE_FLAG = 'basepflag'
#: BOOL (1/2/3); Re-project data from Line of sight, 1 = vertical, 2 = horizontal, 3 = no conversion
REPROJECTION_FLAG = 'prjflag'

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
#: FLOAT; Coordinate in x of reference pixel OR -1 = perform reference pixel search 
REFX = 'refx'
#: FLOAT; Coordinate in y of reference pixel OR -1 = perform reference pixel search 
REFY = 'refy'
#: INT; Number of reference pixel grid search nodes in x dimension
REFNX = "refnx"
#: INT; Number of reference pixel grid search nodes in y dimension
REFNY = "refny"
#: INT; Dimension of reference pixel search window
REF_CHIP_SIZE = 'refchipsize'
#: REAL; Minimum fraction of observations required in reference pixel search window for pixel to be a viable reference pixel
REF_MIN_FRAC = 'refminfrac'

# orbital error correction/parameters
#: BOOL (1/0); Boolean flag controlling whether to apply orbital error correction
ORBITAL_FIT = 'orbfit' 
#: BOOL (1/2); Method for orbital error correction, 1: ifg by ifg/independent, 2: epoch by epoch/network
ORBITAL_FIT_METHOD = 'orbfitmethod'  
#: BOOL (1/2/3) Order of orbital error model, 1 = planar in x and y (2 parameter model, 2 = quadratic in x and y (5 parameter model), 3 = quadratic in x and cubic in y (part-cubic 6 parameter model)
ORBITAL_FIT_DEGREE = 'orbfitdegrees' 
#: INT; Multi look factor for orbital error calculation in x dimension
ORBITAL_FIT_LOOKS_X = 'orbfitlksx'
#: INT; Multi look factor for orbital error calculation in y dimension 
ORBITAL_FIT_LOOKS_Y = 'orbfitlksy'
# ORBITAL_FIT_orbrefest:	 1 BOOLEAN (1/0) # remove reference phase
# ORBITAL_FIT_ orbmaskflag:   1 BOOLEAN (1/0) # mask some patches for orbital correction

# Linear rate/stacking parameters
#: REAL; Threshold ratio between 'model minus observation' residuals and a-priori observation standard deviations for linear rate estimate acceptance (otherwise remove furthest outlier and re-iterate)
LR_NSIG = 'nsig'
#: INT; Number of required input observations per pixel for the linear rate inversion
LR_PTHRESH = 'pthr'
#: REAL; Maximum allowable standard error for pixels in linear rate inversion.
LR_MAXSIG = 'maxsig'

#: INT; Number of required input observations per pixel for time series inversion
TIME_SERIES_PTHRESH = 'ts_pthr'
#: BOOL (1/0); Time series parameter to interpolate across epoch gaps NOT CURRENTLY USED
TIME_SERIES_INTERP = 'ts_interp'


def degree_conv(deg):
    """
    Convenience: convert numerical degree to human readable string
    """

    degree = int(deg)
    if degree == 1:
        return orbital.PLANAR
    if degree == 2:
        return orbital.QUADRATIC
    if degree == 3:
        return orbital.PART_CUBIC
    raise NotImplementedError


def method_conv(meth):
    """
    Convenience: convert numerical method to human readable string
    """

    method = int(meth)
    if method == 1:
        return orbital.INDEPENDENT_METHOD
    if method == 2:
        return orbital.NETWORK_METHOD
    raise NotImplementedError


# Lookup to help convert args to correct type/defaults
# format is	key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = {
    PERP_BASELINE_FLAG : (bool, True),
    AMPLITUDE_FLAG : (bool, False),
    NUM_SETS : (int, 1),
    IFG_CROP_OPT : (int, None), # TODO: default to ALREADY_SAME_SIZE?
    IFG_LKSX : (int, NO_MULTILOOKING),
    IFG_LKSY : (int, NO_MULTILOOKING),
    IFG_XFIRST : (float, None),
    IFG_XLAST : (float, None),
    IFG_YFIRST : (float, None),
    IFG_YLAST : (float, None),
    REPROJECTION_FLAG : (int, 3),

    REFX : (int, -1),
    REFY : (int, -1),
    REFNX : (int, None), # was 50 in original Pirate code
    REFNY : (int, None), # was 50 in original Pirate code
    REF_CHIP_SIZE : (int, None), # defaults to 21 in orig
    REF_MIN_FRAC : (float, 0.8), # uses Pirate default

    ORBITAL_FIT : (bool, True),
    ORBITAL_FIT_METHOD : (method_conv, orbital.NETWORK_METHOD),
    ORBITAL_FIT_DEGREE : (degree_conv, orbital.QUADRATIC),
    ORBITAL_FIT_LOOKS_X : (int, NO_MULTILOOKING),
    ORBITAL_FIT_LOOKS_Y : (int, NO_MULTILOOKING),

    LR_NSIG : (int, 3), # Pirate default
    LR_PTHRESH : (int, 20), # should be based on nepochs since not every project may have 20 epochs
    LR_MAXSIG : (int, 2), # Pirate default

    TIME_SERIES_PTHRESH : (int, None)}


def get_config_params(path):
    """
    Returns a dict for the key:value pairs from the .conf file
    """

    with open(path) as f:
        txt = f.read()

    return _parse_conf_file(txt)


def _parse_conf_file(content):
    """
    Parser for converting text content into a dict of parameters
    """

    def is_valid(line):
        return line != "" and line[0] not in "%#"

    lines = [ln.split() for ln in content.split('\n') if is_valid(ln)]

    # convert "field:   value" lines to [field, value]
    kvpair = [(e[0].rstrip(":"), e[1]) for e in lines if len(e) == 2]
    parameters = dict(kvpair)

    if not parameters:
        raise ConfigException('Cannot parse any parameters from config file')

    return _parse_pars(parameters)


def _parse_pars(pars):
    """
    Parses and converts config file params from text
    """

    for k in PARAM_CONVERSION.keys():
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
    with open(nml) as f:
        return [ln.strip() for ln in f.readlines() if ln != ""]


class ConfigException(Exception):
    """
    Default exception class for configuration errors.
    """

    pass

