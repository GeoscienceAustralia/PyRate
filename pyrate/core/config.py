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
import itertools
import logging
import os
import pathlib
import re
from os.path import splitext, split
# coding: utf-8
from typing import List, Tuple, Dict, Optional

from constants import CONV2TIF, PREPIFG, PROCESS, MERGE
from core.ifgconstants import YEARS_PER_DAY

_logger = logging.getLogger(__name__)

# general constants
NO_MULTILOOKING = 1
ROIPAC = 0
GAMMA = 1
LOG_LEVEL = "INFO"
SIXTEEN_DIGIT_EPOCH_PAIR = r"\d{8}-\d{8}"
TWELVE_DIGIT_EPOCH_PAIR = r"\d{6}-\d{6}"
EIGHT_DIGIT_EPOCH = r"\d{8}"
MINIMUM_NUMBER_EPOCHS = 3

# constants for lookups
#: STR; Name of input interferogram list file
IFG_FILE_LIST = "ifgfilelist"
#: BOOL (0/1); The interferogram processor used (0==ROIPAC, 1==GAMMA)
PROCESSOR = "processor"
#: STR; Name of directory for saving output products
OUT_DIR = "outdir"
#: STR; Name of Digital Elevation Model file
DEM_FILE = "demfile"
#: STR; Name of the DEM header file
DEM_HEADER_FILE = "demHeaderFile"
#: STR; Name of directory containing GAMMA SLC header files
SLC_DIR = "slcFileDir"
#: STR; Name of the file list containing the pool of available SLC headers
SLC_FILE_LIST = "slcfilelist"

# STR; The projection of the input interferograms.
# TODO: only used in tests; deprecate?
INPUT_IFG_PROJECTION = "projection"
#: FLOAT; The no data value in the interferogram files.
NO_DATA_VALUE = "noDataValue"
#: FLOAT; No data averaging threshold for prepifg
NO_DATA_AVERAGING_THRESHOLD = "noDataAveragingThreshold"
# BOOL (1/2/3); Re-project data from Line of sight, 1 = vertical, 2 = horizontal, 3 = no conversion
# REPROJECTION = 'prjflag' # NOT CURRENTLY USED
#: BOOL (0/1): Convert no data values to Nan
NAN_CONVERSION = "nan_conversion"

# Prepifg parameters
#: BOOL (1/2/3/4); Method for cropping interferograms, 1 = minimum overlapping area (intersection), 2 = maximum area (union), 3 = customised area, 4 = all ifgs already same size
IFG_CROP_OPT = "ifgcropopt"
#: INT; Multi look factor for interferogram preparation in x dimension
IFG_LKSX = "ifglksx"
#: INT; Multi look factor for interferogram preparation in y dimension
IFG_LKSY = "ifglksy"
#: FLOAT; Minimum longitude for cropping with method 3
IFG_XFIRST = "ifgxfirst"
#: FLOAT; Maximum longitude for cropping with method 3
IFG_XLAST = "ifgxlast"
#: FLOAT; Minimum latitude for cropping with method 3
IFG_YFIRST = "ifgyfirst"
#: FLOAT; Maximum latitude for cropping with method 3
IFG_YLAST = "ifgylast"

# reference pixel parameters
#: INT; Longitude (decimal degrees) of reference pixel, or if left blank a search will be performed
REFX = "refx"
#: INT; Latitude (decimal degrees) of reference pixel, or if left blank a search will be performed
REFY = "refy"
#: INT; Number of reference pixel grid search nodes in x dimension
REFNX = "refnx"
#: INT; Number of reference pixel grid search nodes in y dimension
REFNY = "refny"
#: INT; Dimension of reference pixel search window (in number of pixels)
REF_CHIP_SIZE = "refchipsize"
#: FLOAT; Minimum fraction of observations required in search window for pixel to be a viable reference pixel
REF_MIN_FRAC = "refminfrac"
#: BOOL (1/2); Reference phase estimation method (1: median of the whole interferogram, 2: median within the window surrounding the reference pixel)
REF_EST_METHOD = "refest"

# coherence masking parameters
#: BOOL (0/1); Perform coherence masking (1: yes, 0: no)
COH_MASK = "cohmask"
#: FLOAT; Coherence threshold for masking
COH_THRESH = "cohthresh"
#: STR; Directory containing coherence files
COH_FILE_DIR = "cohfiledir"
#: STR; Name of the file list containing the pool of available coherence files
COH_FILE_LIST = "cohfilelist"

# atmospheric error correction parameters NOT CURRENTLY USED
APS_CORRECTION = "apscorrect"
APS_METHOD = "apsmethod"
APS_INCIDENCE_MAP = "incidencemap"
APS_INCIDENCE_EXT = "APS_INCIDENCE_EXT"
APS_ELEVATION_MAP = "elevationmap"
APS_ELEVATION_EXT = "APS_ELEVATION_EXT"

# orbital error correction/parameters
#: BOOL (1/0); Perform orbital error correction (1: yes, 0: no)
ORBITAL_FIT = "orbfit"
#: BOOL (1/2); Method for orbital error correction (1: independent, 2: network)
ORBITAL_FIT_METHOD = "orbfitmethod"
#: BOOL (1/2/3) Polynomial order of orbital error model (1: planar in x and y - 2 parameter model, 2: quadratic in x and y - 5 parameter model, 3: quadratic in x and cubic in y - part-cubic 6 parameter model)
ORBITAL_FIT_DEGREE = "orbfitdegrees"
#: INT; Multi look factor for orbital error calculation in x dimension
ORBITAL_FIT_LOOKS_X = "orbfitlksx"
#: INT; Multi look factor for orbital error calculation in y dimension
ORBITAL_FIT_LOOKS_Y = "orbfitlksy"

# Stacking parameters
#: FLOAT; Threshold ratio between 'model minus observation' residuals and a-priori observation standard deviations for stacking estimate acceptance (otherwise remove furthest outlier and re-iterate)
LR_NSIG = "nsig"
#: INT; Number of required observations per pixel for stacking to occur
LR_PTHRESH = "pthr"
#: FLOAT; Maximum allowable standard error for pixels in stacking
LR_MAXSIG = "maxsig"

# atmospheric delay errors fitting parameters NOT CURRENTLY USED
# atmfitmethod = 1: interferogram by interferogram; atmfitmethod = 2, epoch by epoch
# ATM_FIT = 'atmfit'
# ATM_FIT_METHOD = 'atmfitmethod'

# APS correction parameters
#: BOOL (0/1) Perform APS correction (1: yes, 0: no)
APSEST = "apsest"
# temporal low-pass filter parameters
#: INT (1/2/3); Method for temporal filtering (1: Gaussian, 2: Triangular, 3: Mean filter)
TLPF_METHOD = "tlpfmethod"
#: FLOAT; Cutoff time for gaussian filter in years;
TLPF_CUTOFF = "tlpfcutoff"
#: INT; Number of required input observations per pixel for temporal filtering
TLPF_PTHR = "tlpfpthr"
# spatially correlated noise low-pass filter parameters
#: INT (1/2); Method for spatial filtering(1: butterworth; 2: gaussian)
SLPF_METHOD = "slpfmethod"
#: FLOAT; Cutoff  value for both butterworth and gaussian filters in km
SLPF_CUTOFF = "slpfcutoff"
#: INT; Order of butterworth filter (default 1)
SLPF_ORDER = "slpforder"
#: INT (1/0); Do spatial interpolation at NaN locations (1 for interpolation, 0 for zero fill)
SLPF_NANFILL = "slpnanfill"
#: #: STR; Method for spatial interpolation (one of: linear, nearest, cubic), only used when slpnanfill=1
SLPF_NANFILL_METHOD = "slpnanfill_method"

# Time series parameters
#: BOOL (1/0); Perform time series calculation (1: yes, 0: no)
TIME_SERIES_CAL = "tscal"
#: INT (1/2); Method for time series inversion (1: Laplacian Smoothing; 2: SVD)
TIME_SERIES_METHOD = "tsmethod"
#: INT; Number of required input observations per pixel for time series inversion
TIME_SERIES_PTHRESH = "ts_pthr"
#: INT (1/2); Order of Laplacian smoothing operator, first or second order
TIME_SERIES_SM_ORDER = "smorder"
#: FLOAT; Laplacian smoothing factor (values used is 10**smfactor)
TIME_SERIES_SM_FACTOR = "smfactor"
# tsinterp is automatically assigned in the code; not needed in conf file
# TIME_SERIES_INTERP = 'tsinterp'

#: BOOL (0/1/2); Use parallelisation/Multi-threading (0: in serial, 1: in parallel by rows, 2: in parallel by pixel)
PARALLEL = "parallel"
#: INT; Number of processes for multi-threading
PROCESSES = "processes"

# Orbital error correction constants for conversion to readable strings
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2
PLANAR = 1
QUADRATIC = 2
PART_CUBIC = 3

# Orbital error name look up for logging
ORB_METHOD_NAMES = {INDEPENDENT_METHOD: "INDEPENDENT", NETWORK_METHOD: "NETWORK"}
ORB_DEGREE_NAMES = {PLANAR: "PLANAR", QUADRATIC: "QUADRATIC", PART_CUBIC: "PART CUBIC"}

# dir for temp files
TMPDIR = "tmpdir"

# CONFIG UTILS - TO BE MOVED?
def parse_namelist(nml):
    """Parses name list file into array of paths

    Args:
      str: nml: interferogram file list
      nml: returns: list of interferogram file names

    Returns:
      list: list of interferogram file names

    """
    with open(nml) as f_in:
        lines = [line.rstrip() for line in f_in]
    return filter(None, lines)


def write_config_file(params, output_conf_file):
    """Takes a param object and writes the config file. Reverse of get_conf_params.
    
    Args:
      dict: params: parameter dictionary
      str: output_conf_file: output file name

    Args:
      output_conf_file: 
      params: 

    Returns:
      list: config file

    """
    with open(output_conf_file, "w") as f:
        for k, v in params.items():
            if v is not None:
                f.write("".join([k, ":\t", str(v), "\n"]))
            else:
                f.write("".join([k, ":\t", "", "\n"]))


def coherence_paths_for(path, params, tif=False) -> str:
    """Returns path to coherence file for given interferogram. Pattern matches
    based on epoch in filename.
    
    Example:
        '20151025-20160501_eqa_filt.cc'
        Date pair is the epoch.
    
    Args:
        path: Path to intergerogram to find coherence file for.
    
    Args:
      tif: Find converted tif if True (Default value = False)
      path: param params:

    Args:
      path: 
      params: 
      tif:  (Default value = False)

    Returns:
      Path to coherence file.

    """
    _, filename = split(path)
    pattern = re.compile(r"\d{8}-\d{8}")
    epoch = re.match(pattern, filename).group(0)
    coherence_dir = params[COH_FILE_DIR]
    matches = [name for name in parse_namelist(params[COH_FILE_LIST]) if epoch in name]
    if tif:
        names_exts = [os.path.splitext(m) for m in matches]
        matches = [ne[0] + ne[1].replace(".", "_") + ".tif" for ne in names_exts]

    for i, m in enumerate(matches):
        if "_tif" in m:
            matches[i] = m.replace("_tif", "")

    return [os.path.join(coherence_dir, m) for m in matches]


def coherence_paths(params) -> List[str]:
    """Returns paths to corresponding coherence files for given IFGs. Assumes
    that each IFG has a corresponding coherence file in the coherence file
    directory and they share epoch prefixes.
    
    Args:
      ifg_paths: List of paths to intergerogram files.

    Args:
      params: 

    Returns:
      A list of full paths to coherence files.

    """
    ifg_file_list = params.get(IFG_FILE_LIST)
    ifgs = parse_namelist(ifg_file_list)
    paths = [coherence_paths_for(ifg, params) for ifg in ifgs]
    return list(itertools.chain.from_iterable(paths))


def mlooked_path(path, looks, crop_out):
    """Adds suffix to ifg path, for creating a new path for multilooked files.

    Args:
      str: path: original interferogram path
      int: looks: number of range looks applied
      int: crop_out: crop option applied
      path: param looks:
      crop_out: returns: multilooked file name
      looks: 

    Returns:
      str: multilooked file name

    """
    base, ext = splitext(path)
    return "{base}_{looks}rlks_{crop_out}cr{ext}".format(base=base, looks=looks, crop_out=crop_out, ext=ext)

class ConfigException(Exception):
    """Default exception class for configuration errors."""
