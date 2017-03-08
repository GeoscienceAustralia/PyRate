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
This Python module implements an algorithm to search for the location
of the interferometric reference pixel
"""
import os
import logging
from itertools import product
import numpy as np
from numpy import isnan, std, mean, sum as nsum
from joblib import Parallel, delayed

import pyrate.config as cf
from pyrate.shared import Ifg

log = logging.getLogger(__name__)


# TODO: move error checking to config step (for fail fast)
def ref_pixel(ifgs, params):
    """
    Returns (y,x) reference pixel coordinate from given ifgs.

    If the config file REFX or REFY values are empty or subzero, the search for
    the reference pixel is performed. If the REFX|Y values are within the
    bounds of the raster, a search is not performed. REFX|Y values outside
    the upper bounds cause an exception.

    Parameters
    ----------
    ifgs: list
        sequence of interferograms.
    params: dict
        parameters dictionary

    """
    half_patch_size, thresh, grid = ref_pixel_setup(ifgs, params)
    parallel = params[cf.PARALLEL]
    if parallel:
        phase_data = [i.phase_data for i in ifgs]
        mean_sds = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(ref_pixel_multi)(g, half_patch_size, phase_data,
                                     thresh, params) for g in grid)
        refy, refx = filter_means(mean_sds, grid)
    else:
        phase_data = [i.phase_data for i in ifgs]
        mean_sds = []
        for g in grid:
            mean_sds.append(ref_pixel_multi(
                g, half_patch_size, phase_data, thresh, params))
        refy, refx = filter_means(mean_sds, grid)

    if refy and refx:
        return refy, refx

    raise RefPixelError("Could not find a reference pixel")


def filter_means(mean_sds, grid):
    """
    filter the mean standard deviations from each ref pixel block
    Parameters
    ----------
    mean_sds: list
        list of mean standard deviations from each ref pixel grid
    grid: list
        list of grid tuples
    Return
    ------
    tuple of (refy, refx)
    """
    log.info('Filtering means during reference pixel computation')
    refp_index = np.nanargmin(mean_sds)
    return grid[refp_index]


def ref_pixel_setup(ifgs_or_paths, params):
    """
    sets up the grid for reference pixel computation
    Also saves numpy files for later use during ref pixel computation
    """
    log.info('Setting up ref pixel computation')
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
    validate_chipsize(chipsize, head)
    validate_minimum_fraction(min_frac)
    validate_search_win(refnx, refny, chipsize, head)
    # pre-calculate useful amounts
    half_patch_size = chipsize // 2
    chipsize = half_patch_size * 2 + 1
    thresh = min_frac * chipsize * chipsize
    # do window searches across dataset, central pixel of stack with smallest
    # mean is the reference pixel
    rows, cols = head.shape
    ysteps = step(rows, refny, half_patch_size)
    xsteps = step(cols, refnx, half_patch_size)
    log.info('Ref pixel setup finished')
    return half_patch_size, thresh, list(product(ysteps, xsteps))


def ref_pixel_mpi(process_grid, half_patch_size, ifgs, thresh, params):
    """
    convenience function for ref pixel calculation for each process
    """
    log.info('Ref pixel calculation started')
    mean_sds = []
    for g in process_grid:
        mean_sds.append(ref_pixel_multi(g, half_patch_size, ifgs, thresh,
                                        params))
    return mean_sds


def ref_pixel_multi(g, half_patch_size, phase_data_or_ifg_paths,
                    thresh, params):
    """
    convenience function for for ref pixel optimisation
    """
    # pylint: disable=invalid-name
    # phase_data_or_ifg is list of ifgs
    y, x, = g
    if isinstance(phase_data_or_ifg_paths[0], str):
        # this consumes a lot less memory
        # one ifg.phase_data in memory at any time
        data = []
        output_dir = params[cf.OUT_DIR]
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


def step(dim, ref, radius):
    """
    Helper func: returns xrange obj of axis indicies for a search window.

    :param dim: total length of the grid dimension.
    :param ref: the desired number of steps.
    :param radius: the number of cells from the centre of the chip eg.
    (chipsize / 2).
    """

    # if ref == 1:
    #     # centre a single search step
    #     return xrange(dim // 2, dim, dim)  # fake step to ensure single xrange value

    # if ref == 2: # handle 2 search windows, method below doesn't cover the case
    #     return [radius, dim-radius-1]
    # max_dim = dim - (2*radius)  # max possible number for refn(x|y)
    # step = max_dim // (ref-1)
    step_size = dim // ref  # same as in Matlab
    return range(radius, dim-radius, step_size)


def validate_chipsize(chipsize, head):
    """sanity check min chipsize"""
    if chipsize is None:
        raise cf.ConfigException('Chipsize is None')

    if chipsize < 3 or chipsize > head.ncols or (chipsize % 2 == 0):
        msg = "Chipsize setting must be >=3 and at least <= grid width"
        raise ValueError(msg)
    log.info('Chipsize validation successful')


def validate_minimum_fraction(min_frac):
    """sanity check min fraction"""
    if min_frac is None:
        raise cf.ConfigException('Minimum fraction is None')

    if min_frac < 0.0 or min_frac > 1.0:
        raise ValueError("Minimum fraction setting must be >= 0.0 and <= 1.0 ")


def validate_search_win(refnx, refny, chipsize, head):
    """sanity check X|Y steps"""
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
