import re
import itertools
import logging

from pyrate.core.ifgconstants import YEARS_PER_DAY
from pyrate.core import config
from pyrate.core.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR,
    SLC_FILE_LIST,
    COH_MASK,
    COH_THRESH,
    COH_FILE_DIR,
    COH_FILE_LIST,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    REFNX,
    REFNY,
    REF_CHIP_SIZE,
    REF_MIN_FRAC,
    ORBITAL_FIT,
    ORBITAL_FIT_METHOD,
    ORBITAL_FIT_DEGREE,
    ORBITAL_FIT_LOOKS_X,
    ORBITAL_FIT_LOOKS_Y,
    LR_NSIG,
    LR_MAXSIG,
    LR_PTHRESH,
    APSEST,
    TLPF_METHOD,
    TLPF_CUTOFF,
    TLPF_PTHR,
    SLPF_METHOD,
    SLPF_CUTOFF,
    SLPF_ORDER,
    SLPF_NANFILL,
    TIME_SERIES_CAL,
    TIME_SERIES_PTHRESH,
    TIME_SERIES_SM_FACTOR,
    TIME_SERIES_SM_ORDER,
    TIME_SERIES_METHOD,
    PARALLEL,
    PROCESSES,
    NAN_CONVERSION,
    NO_DATA_AVERAGING_THRESHOLD,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
    APS_METHOD,
    APS_CORRECTION,
    ConfigException)

16_DIGIT_EPOCH_PAIR_PATTERN = r'\d{8}-\d{8}'
12_DIGIT_EPOCH_PAIR_PATTERN = r'\d(6)-\d{6}'
8_DIGIT_EPOCH_PATTERN = r'\d{8}'


PARAM_VALIDATION = {
    OBS_DIR: (
        lambda a: a is not None and os.path.exists(a),
        f"'{OBS_DIR}': directory must be provided and must exist."
    ),
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
        lambda a: a == 1 or a == 2 or a == 3 or a == 4,
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
        lambda a: a == 0 or a == 1,
        f"'{COH_MASK}': must select option 0 or 1."
    ),
    REFX: (
        lambda a: True,
        "Any int value valid."
    ),
    REFY: (
        lambda a: True, 
        "Any int value valid."
    ),
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
        lambda a: a == 1 or a == 2,
        f"'{REF_EST_METHOD}': must select option 1 or 2."
    ), 
    ORBITAL_FIT: (
        lambda a: a == 0 or a == 1, 
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
        lambda a: a == 0 or a == 1,
        f"'{APSEST}': must select option 0 or 1." 
    ),
    TIME_SERIES_CAL: (
        lambda a: a == 0 or a == 1, 
        f"'{TIME_SERIES_CAL}': must select option 0 or 1."
    ),
    PARALLEL: (
        lambda a: a == 0 or a == 1 or a == 2, 
        f"'{PARALLEL}': must select option 0 or 1 or 2."
    ),
    PROCESSES: (
        lambda a: a >= 1,
        f"'{PROCESSES}': must be >= 1."
    ),
    PROCESSOR: (
        lambda a: a == 0 or a == 1,
        f"'{PROCESSOR}': must select option 0 or 1."
    ),
    NAN_CONVERSION: (
        lambda a: a == 0 or a == 1, 
        f"'{NAN_CONVERSION}': must select option 0 or 1."
    ),
    NO_DATA_AVERAGING_THRESHOLD: (
        lambda a: True, 
        "Any float value valid."),
}
"""dict: basic validation functions for compulsory parameters."""

CUSTOM_CROP_VALIDATION = {
    IFG_XFIRST: (
        lambda a: is not None, 
        f"'{IFG_XFIRST}': must be provided."
    ),
    IFG_XLAST: (
        lambda a: is not None, 
        f"'{IFG_XLAST}': must be provided."
    ),
    IFG_YFIRST: (
        lambda a: is not None,
        f"'{IFG_YFIRST}': must be provided."
    ),
    IFG_YLAST: (
        lambda a: is not None,
        f"'{IFG_YLAST}': must be provided.."
    ),
}
"""dict: basic validation functions for custom cropping parameters."""

GAMMA_VALIDATION = {
    SLC_DIR: (
        lambda a: os.path.exists(a) if a is not None else True,
        f"'{SLC_DIR}': directory must must exist."
    ),
    SLC_FILE_LIST: (
        lambda a: a is not None and os.path.exists(a),
        f"'{SLC_FILE_LIST}': file must be provided and must exist."
    ),
}
"""dict: basic validation functions for gamma parameters."""

COHERENCE_VALIDATION = {
    COH_THRESH: (
        lambda a: 0.0 <= a <= 1.0,
        f"'{COH_THRESH}': must be between 0.0 and 1.0 (inclusive)."
    ),
    COH_FILE_DIR: (
        lambda a: os.path.exists(a) if a is not None else True,
        f"'{COH_FILE_DIR}': directory must exist."
    ),
    COH_FILE_LIST: (
        lambda a: a is not None and os.path.exists(a),
        f"'{COH_FILE_LIST}': file must be provided and must exist."
    ),
}
"""dict: basic validation functions for coherence parameters."""

ORBITAL_FIT_VALIDATION = {
    ORBITAL_FIT_METHOD: (
        lambda a: a == 1 or a == 2 , 
        f"'{ORBITAL_FIT_METHOD}': must select option 1 or 2."
    ),
    ORBITAL_FIT_DEGREE: (
        lambda a: a == 1 or a == 2 or a == 3, 
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

APSEST_VALIDATION = {
    TLPF_METHOD: (
        lambda a: a == 1 or a == 2 or a == 3, 
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
        lambda a: a == 1 or a == 2,
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
        lambda a: a == 0 or a == 1, 
        f"'{SLPF_NANFILL}': must select option 0 or 1."
    ),
}
"""dict: basic validation functions for atmospheric correction parameters."""

TIME_SERIES_VALIDATION = {
    TIME_SERIES_PTHRESH: (
        lambda a: a >= 1, 
        f"'{TIME_SERIES_PTHRESH}': must be >= 1"
    ),
    #TODO: Matt to investigate smoothing factor values.
    TIME_SERIES_SM_FACTOR: (
        lambda a: True, 
        f"'{TIME_SERIES_SM_FACTOR}':"
    ),
    TIME_SERIES_SM_ORDER: (
        lambda a: a == 1 or a == 2,
         f"'{TIME_SERIES_SM_ORDER}': must select option 1 or 2." 
    ),
    TIME_SERIES_METHOD: (
        lambda a: a == 1 or a == 2,
        f"'{TIME_SERIES_METHOD}': must select option 1 or 2."
    ),
}
"""dict: basic vaidation functions for time series parameters."""

def validate_parameters(pars):
    """
    Calls validation functions on each parameter, collects errors and 
    raises an exception with collected errors if errors occurred.

    Args:
        pars: the parameters dictionary.

    Raises:
        ConfigException: if errors occur during parameter validation.
    """
    GAMMA = pars[PROCESSOR] == 1
    ifl = pars[IFG_FILE_LIST]

    validate_compulsory_parameters(pars)
    validate_optional_parameters(pars)

    if GAMMA:
        validate_epochs(ifl, 16_DIGIT_EPOCH_PAIR_PATTERN)
    else:
        validate_epochs(ifl, 12_DIGIT_EPOCH_PAIR_PATTERN)

    #validate_minimum_epochs()
    validate_ifgs(ifl, pars[OBS_DIR])
    validate_obs_thresholds(ifl, pars)
    #validate_epoch_thresholds()

    if GAMMA:
        validate_epochs(pars[SLC_FILE_LIST], 8_DIGIT_EPOCH_PATTERN)
        validate_gamma_headers(ifl, pars[SLC_FILE_LIST], pars[SLC_DIR])

    if pars[COH_MASK]:
        validate_epochs(pars[COH_FILE_LIST], 16_DIGIT_EPOCH_PATTERN)
        validate_coherence_files(ifl, pars)

def _raise_errors(errors):
    if errors:
        errors.insert(0, "invalid parameters")
        raise ConfigException('\n'.join(errors))
    else:
        return True

def validate_compulsory_parameters(pars):
    # Basic validation of parameters that are always used.
    errors = []
    for k in pars.keys():
        validator = PARAM_VALIDATION.get(k)
        if validator is None:
            _logger.debug(f"No validator implemented for '{k}'.")
            continue
        if not validator[0](pars[k]):
            errors.append(validator[1])

    return _raise_errors(errors)
    
def validate_optional_parameters(pars):
    def validate(on, key, validators, pars):
        errors = []
        if on:
            for k, validator in validators.items():
                if not validator[0](pars[k]):
                    errors.append(validator[1])
        return errors

    errors = []

    errors.extend(
        validate(COH_MASK, COHERENCE_VALIDATION, pars))
    errors.extend(
        validate(APSEST, APSEST_VALIDATION, pars))
    errors.extend(
        validate(TIME_SERIES_CAL, TIME_SERIES_VALIDATION, pars))
    errors.extend(
        validate(ORBITAL_FIT, ORBITAL_FIT_VALIDATION, pars))
    errors.extend(
        validate(pars[PROCESSOR]==GAMMA, GAMMA_VALIDATION, pars))
    errors.extent(
        validate(pars[IFG_CROP_OPT]==3, CUSTOM_CROP_VALIDATION, pars))

    return _raise_errors(errors)
     
def validate_epochs(file_list, pattern):
    errors = []
    PTN = re.compile(pattern)
    filenames = config.parse_namelist(file_list)
    for fn in filenames:
        epochs = PTN.findall(fn)
        if len(epochs) == 0:
            errors.append(f"'{file_list}': {fn} does not contain an epoch of "
                           "format {pattern}."
        if len(epochs) > 1:
            errors.append(f"'{file_list}': {fn} does contains more than "
                           "one epoch of {pattern}. There must be only one "
                           "epoch in the filename.")

    return _raise_errors(errors)

def validate_ifgs(ifg_file_list, obs_dir):
    """
    Checks that the IFGs specified in IFG_FILE_LIST exist.

    Args:
        ifgs: list of IFG filenames.
        obs_dir: the observations directory where the IFGs should exist.

    Returns:
        A list of error messages with one for each IFG that doesn't exist, 
        otherwise an empty list if all IFGs exist.
    """
    errors = []
    ifgs = config.parse_namelist(ifg_file_list)
    ifg_paths = [os.path.join(obs_dir, ifg) for ifg in ifgs]
    for path in ifg_paths:
        if not os.path.exists(path):
            fname = os.path.split(path)[1]
            errors.append(f"'{IFG_FILE_LIST}': interferogram '{fname}' does not exist.")

    return _raise_errors(errors)

def validate_obs_thresholds(ifg_file_list, pars):
    """
    Validates parameters that specify an observations threshold.

    Args:
        n_ifgs: the number of IFGs specified in IFG_FILE_LIST.
        pars: the parameters dictionary.
        key: key for the observations threshold being validated.

    Returns:
        An error message if n_ifgs is less than the value set by the 
        threshold parameterm, otherwise None.
    """
    def validate(n, p, k):
        thresh = p[k]
        if thresh > n:
            return [f"'{key}': not enough interferograms have been specified "
                    f"({n_ifgs}) to satisfy threshold ({thresh})."]
        else:
            return []
    
    errors = []
    n_ifgs = len(list(config.parse_namelist(ifg_file_list)))
    errors.extend(validate(n_ifgs, pars, TIME_SERIES_PTHRESH))
    if pars[APSEST]:
        errors.extend(validate(n_ifgs, pars, TLPF_PTHR))

    return _raise_errors(errors)
            
def validate_gamma_headers(ifg_file_list, slc_file_list, slc_dir):
    """
    Validates that gamme headers exist for specified IFGs.

    Args:
        ifgs: list of IFG filenames.
        slc_dir: the slc directory where gamma headers should exist.

    Returns:
        A list of error messages with one for each IFG that doesn't have 
        2 matching headers (one for each epoch), otherwise an empty list if
        all headers exist.
    """
    from pyrate.core.gamma import get_header_paths
    errors = []
    
    for ifg in config.parse_namelist(ifg_file_list):
        headers = get_header_paths(ifg, slc_file_list, slc_dir)
        if len(headers) < 2:
            errors.append(f"'{SLC_DIR}': Headers not found for interferogram "
                           "'{ifg}'.")

    return _raise_errrors(errors)

def validate_coherence_files(ifg_file_list, pars):
    errors = []
    
    for ifg in config.parse_namelist(ifg_file_list):
        paths = coherence_paths_for(ifg, pars)
        if len(paths) == 0:
            errors.append(f"'{COH_FILE_DIR}': no coherence files found for "
                          f"intergerogram '{ifg}'.")
        elif len(paths) > 2:
            errors.append(f"'{COH_FILE_DIR}': found more than one coherence "
                          f"file for '{ifg}'. There must be only one "
                          f"coherence file per interferogram. Found {paths}.")

    return _raise_errors(errors)

