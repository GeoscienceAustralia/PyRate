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
This Python module implements an algorithm to search for the location
of the interferometric reference pixel
"""
import os
from os.path import join
from itertools import product

import numpy as np
from numpy import isnan, std, mean, sum as nsum
from joblib import Parallel, delayed

import pyrate.core.config as cf
from pyrate.core.shared import Ifg
from pyrate.core.shared import joblib_log_level
from pyrate.core.logger import pyratelogger as log


# TODO: move error checking to config step (for fail fast)
# TODO: this function is not used. Plan removal
def ref_pixel(ifgs, params):
    """
    Determines the most appropriate reference pixel coordinate by conducting
    a grid search and calculating the mean standard deviation with patches
    around candidate pixels from the given interferograms.

    If the config file REFX or REFY values are empty or negative, the search
    for the reference pixel is performed. If the REFX|Y values are within the
    bounds of the raster, a search is not performed. REFX|Y values outside
    the upper bounds cause an exception.

    :param list ifgs: List of interferogram objects
    :param dict params: Dictionary of configuration parameters

    :return: tuple of best REFX and REFY coordinates
    :rtype: tuple
    """
    half_patch_size, thresh, grid = ref_pixel_setup(ifgs, params)
    parallel = params[cf.PARALLEL]
    if parallel:
        phase_data = [i.phase_data for i in ifgs]
        mean_sds = Parallel(n_jobs=params[cf.PROCESSES], 
                            verbose=joblib_log_level(cf.LOG_LEVEL))(
            delayed(_ref_pixel_multi)(g, half_patch_size, phase_data,
                                     thresh, params) for g in grid)
        refxy = find_min_mean(mean_sds, grid)
    else:
        phase_data = [i.phase_data for i in ifgs]
        mean_sds = []
        for g in grid:
            mean_sds.append(_ref_pixel_multi(g, half_patch_size, phase_data, thresh, params))
        refxy = find_min_mean(mean_sds, grid)

    if isinstance(refxy, ValueError):
        raise ValueError('Refpixel calculation not possible!')

    refy, refx = refxy

    if refy and refx:
        return refy, refx

    raise RefPixelError("Could not find a reference pixel")


def find_min_mean(mean_sds, grid):
    """
    Determine the ref pixel block with minimum mean value

    :param list mean_sds: List of mean standard deviations from each
        reference pixel grid
    :param list grid: List of ref pixel coordinates tuples

    :return: Tuple of (refy, refx) with minimum mean
    :rtype: tuple    
    """
    log.debug('Ranking ref pixel candidates based on mean values')
    try:
        refp_index = np.nanargmin(mean_sds)
        return grid[refp_index]
    except ValueError as v:
        log.error(v)
        return v


def ref_pixel_setup(ifgs_or_paths, params):
    """
    Sets up the grid for reference pixel computation and saves numpy files
    to disk for later use during ref pixel computation.
        
    :param list ifgs_or_paths: List of interferogram filenames or Ifg objects
    :param dict params: Dictionary of configuration parameters
    
    :return: half_patch_size: size of patch
    :rtype: float
    :return: thresh
    :rtype: float
    :return: list(product(ysteps, xsteps))
    :rtype: list
    """
    log.debug('Setting up ref pixel computation')
    refnx, refny, chipsize, min_frac = params[cf.REFNX], \
                                       params[cf.REFNY], \
                                       params[cf.REF_CHIP_SIZE], \
                                       params[cf.REF_MIN_FRAC]
    if len(ifgs_or_paths) < 1:
        msg = 'Reference pixel search requires 2+ interferograms'
        raise RefPixelError(msg)

    if isinstance(ifgs_or_paths[0], str):
        head = Ifg(ifgs_or_paths[0])
        head.open(readonly=True)
    else:
        head = ifgs_or_paths[0]

    # sanity check inputs
    _validate_chipsize(chipsize, head)
    _validate_minimum_fraction(min_frac)
    _validate_search_win(refnx, refny, chipsize, head)
    # pre-calculate useful amounts
    half_patch_size = chipsize // 2
    chipsize = half_patch_size * 2 + 1
    thresh = min_frac * chipsize * chipsize
    # do window searches across dataset, central pixel of stack with smallest
    # mean is the reference pixel
    rows, cols = head.shape
    ysteps = _step(rows, refny, half_patch_size)
    xsteps = _step(cols, refnx, half_patch_size)
    log.debug('Ref pixel setup finished')
    return half_patch_size, thresh, list(product(ysteps, xsteps))


def save_ref_pixel_blocks(grid, half_patch_size, ifg_paths, params):
    """
    Save reference pixel grid blocks to numpy array files on disk

    :param list grid: List of tuples (y, x) corresponding to ref pixel grids
    :param int half_patch_size: patch size in pixels
    :param list ifg_paths: list of interferogram paths
    :param dict params: Dictionary of configuration parameters

    :return: None, file saved to disk
    """
    log.debug('Saving ref pixel blocks')
    outdir = params[cf.TMPDIR]
    for pth in ifg_paths:
        ifg = Ifg(pth)
        ifg.open(readonly=True)
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.convert_to_mm()
        for y, x in grid:
            data = ifg.phase_data[y - half_patch_size:y + half_patch_size + 1,
                                  x - half_patch_size:x + half_patch_size + 1]

            data_file = join(outdir, 'ref_phase_data_{b}_{y}_{x}.npy'.format(
                    b=os.path.basename(pth).split('.')[0], y=y, x=x))
            np.save(file=data_file, arr=data)
        ifg.close()
    log.debug('Saved ref pixel blocks')


def _ref_pixel_mpi(process_grid, half_patch_size, ifgs, thresh, params):
    """
    Convenience function for MPI-enabled ref pixel calculation
    """
    log.debug('Ref pixel calculation started')
    mean_sds = []
    for g in process_grid:
        mean_sds.append(_ref_pixel_multi(g, half_patch_size, ifgs, thresh, params))
    return mean_sds


def _ref_pixel_multi(g, half_patch_size, phase_data_or_ifg_paths,
                    thresh, params):
    """
    Convenience function for ref pixel optimisation
    """
    # pylint: disable=invalid-name
    # phase_data_or_ifg is list of ifgs
    y, x, = g
    if isinstance(phase_data_or_ifg_paths[0], str):
        # this consumes a lot less memory
        # one ifg.phase_data in memory at any time
        data = []
        output_dir = params[cf.TMPDIR]
        for p in phase_data_or_ifg_paths:
            data_file = os.path.join(output_dir,
                                     'ref_phase_data_{b}_{y}_{x}.npy'.format(
                                         b=os.path.basename(p).split('.')[0],
                                         y=y, x=x))
            data.append(np.load(file=data_file))
    else:  # phase_data_or_ifg is phase_data list
        data = [p[y - half_patch_size:y + half_patch_size + 1,
                  x - half_patch_size:x + half_patch_size + 1]
                for p in phase_data_or_ifg_paths]
    valid = [nsum(~isnan(d)) > thresh for d in data]
    if all(valid):  # ignore if 1+ ifgs have too many incoherent cells
        sd = [std(i[~isnan(i)]) for i in data]
        return mean(sd)
    else:
        return np.nan


def _step(dim, ref, radius):
    """
    Helper: returns range object of axis indices for a search window.

    :param int dim: Total length of the grid dimension
    :param int ref: The desired number of steps
    :param float radius: The number of cells from the centre of the chip eg.
        (chipsize / 2)

    :return: range object of axis indices
    :rtype: range
    """

    # if ref == 1:
    #     # centre a single search step
    #     return xrange(dim // 2, dim, dim)  # fake step to ensure single xrange value

    # if ref == 2: # handle 2 search windows, method below doesn't cover the case
    #     return [radius, dim-radius-1]
    # max_dim = dim - (2*radius)  # max possible number for refn(x|y)
    # step = max_dim // (ref-1)
    step_size = dim // ref
    return range(radius, dim-radius, step_size)


def _validate_chipsize(chipsize, head):
    """
    Sanity check min chipsize
    """
    if chipsize is None:
        raise cf.ConfigException('Chipsize is None')

    if chipsize < 3 or chipsize > head.ncols or (chipsize % 2 == 0):
        msg = "Chipsize setting must be >=3 and at least <= grid width"
        raise ValueError(msg)
    log.debug('Chipsize validation successful')


def _validate_minimum_fraction(min_frac):
    """
    Sanity check min fraction
    """
    if min_frac is None:
        raise cf.ConfigException('Minimum fraction is None')

    if min_frac < 0.0 or min_frac > 1.0:
        raise ValueError("Minimum fraction setting must be >= 0.0 and <= 1.0 ")


def _validate_search_win(refnx, refny, chipsize, head):
    """
    Sanity check X|Y steps
    """
    if refnx is None:
        raise cf.ConfigException('refnx is None')

    max_width = (head.ncols - (chipsize-1))
    if refnx < 1 or refnx > max_width:
        msg = "Invalid refnx setting, must be > 0 and <= %s"
        raise ValueError(msg % max_width)

    if refny is None:
        raise cf.ConfigException('refny is None')

    max_rows = (head.nrows - (chipsize-1))
    if refny < 1 or refny > max_rows:
        msg = "Invalid refny setting, must be > 0 and <= %s"
        raise ValueError(msg % max_rows)


class RefPixelError(Exception):
    '''
    Generic exception for reference pixel errors.
    '''
