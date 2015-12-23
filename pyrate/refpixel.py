'''
Functions for finding the reference pixel in PyRate.

.. codeauthor:: Ben Davies
'''

import pyrate.config as config
from numpy import array, isnan, std, mean, sum as nsum
import numpy as np


# TODO: move error checking to config step (for fail fast)
def ref_pixel(ifgs, refnx, refny, chipsize, min_frac):
    """
    Returns (y,x) reference pixel coordinate from given ifgs.

    If the config file REFX or REFY values are empty or subzero, the search for
    the reference pixel is performed. If the REFX|Y values are within the bounds
    of the raster, a search is not performed. REFX|Y values outside the upper
    bounds cause an exception.

    :param ifgs: sequence of interferograms.
    """

    if len(ifgs) < 1:
        msg = 'Reference pixel search requires 2+ interferograms'
        raise RefPixelError(msg)

    # sanity check inputs
    head = ifgs[0]
    validate_chipsize(chipsize, head)
    validate_minimum_fraction(min_frac)
    validate_search_win(refnx, refny, chipsize, head)

    # pre-calculate useful amounts
    half_patch_size = chipsize // 2
    chipsize = half_patch_size * 2 + 1
    phase_stack = array([i.phase_data for i in ifgs])  # TODO: mem efficiencies?
    thresh = min_frac * chipsize * chipsize
    min_sd = np.finfo(np.float64).max

    # do window searches across dataset, central pixel of stack with smallest
    # mean is the reference pixel
    refx = None
    refy = None

    rows, cols = ifgs[0].shape

    for y in _step(rows, refny, half_patch_size):
        for x in _step(cols, refnx, half_patch_size):
            data = phase_stack[:, y-half_patch_size:y+half_patch_size+1,
                   x-half_patch_size:x+half_patch_size+1]
            valid = [nsum(~isnan(i)) > thresh for i in data]

            if all(valid):  # ignore if 1+ ifgs have too many incoherent cells
                sd = [std(i[~isnan(i)]) for i in data]
                mean_sd = mean(sd)
                if mean_sd < min_sd:
                    min_sd = mean_sd
                    refy, refx = y, x

    if refy and refx:
        return refy, refx

    raise RefPixelError("Could not find a reference pixel")


def _step(dim, ref, radius):
    '''
    Helper func: returns xrange obj of axis indicies for a search window.

    :param dim: total length of the grid dimension.
    :param ref: the desired number of steps.
    :param radius: the number of cells from the centre of the chip eg. (chipsize / 2).
    '''

    # if ref == 1:
    #     # centre a single search step
    #     return xrange(dim // 2, dim, dim)  # fake step to ensure single xrange value

    # if ref == 2: # handle 2 search windows, method below doesn't cover the case
    #     return [radius, dim-radius-1]
    # max_dim = dim - (2*radius)  # max possible number for refn(x|y)
    # step = max_dim // (ref-1)
    step = dim // ref  # same as in Matlab
    return xrange(radius, dim-radius, step)


def validate_chipsize(chipsize, head):
    if chipsize is None:
        raise config.ConfigException('Chipsize is None')

    if chipsize < 3 or chipsize > head.ncols or (chipsize % 2 == 0):
        msg = "Chipsize setting must be >=3 and at least <= grid width"
        raise ValueError(msg)


def validate_minimum_fraction(min_frac):
    if min_frac is None:
        raise config.ConfigException('Minimum fraction is None')

    if min_frac < 0.0 or min_frac > 1.0:
        raise ValueError("Minimum fraction setting must be >= 0.0 and <= 1.0 ")


def validate_search_win(refnx, refny, chipsize, head):
    # sanity check X|Y steps
    if refnx is None:
        raise config.ConfigException('refnx is None')

    max_width = (head.ncols - (chipsize-1))
    if refnx < 1 or refnx > max_width:
        msg = "Invalid refnx setting, must be > 0 and <= %s"
        raise ValueError(msg % max_width)

    if refny is None:
        raise config.ConfigException('refny is None')

    max_rows = (head.nrows - (chipsize-1))
    if refny < 1 or refny > max_rows:
        msg = "Invalid refny setting, must be > 0 and <= %s"
        raise ValueError(msg % max_rows)


class RefPixelError(Exception):
    '''
    Generic exception for reference pixel errors.
    '''

    pass
