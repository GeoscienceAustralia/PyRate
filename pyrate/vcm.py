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
This Python module implements covariance calculation and
Variance/Covariance matrix functionality. The algorithms
are based on Matlab codes 'cvdcalc.m' and 'vcmt.m' from
the Pirate package.
"""
from __future__ import print_function
from numpy import array, where, isnan, real, imag, sqrt, meshgrid
from numpy import zeros, vstack, ceil, mean, exp, reshape
from numpy.linalg import norm
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.optimize import fmin

from pyrate import shared
from pyrate.shared import PrereadIfg
from pyrate.algorithm import master_slave_ids


def pendiffexp(alphamod, cvdav):
    """
    Fits an exponential model to data.

    :param float alphamod: Exponential decay exponent.
    :param array cvdav: Function magnitude at 0 radius (2 col array of radius,
    variance)
    """
    # pylint: disable=invalid-name

    # maxvar usually at zero lag
    mx = cvdav[1, 0]
    return norm(cvdav[1, :] - (mx * exp(-alphamod * cvdav[0, :])))


# this is not used any more
def unique_points(points):
    """
    Returns unique points from a list of coordinates.

    :param points: Sequence of (y,x) or (x,y) tuples.
    """
    return vstack([array(u) for u in set(points)])


def cvd(ifg_path, params, calc_alpha=False):
    """
    Calculate average covariance versus distance (autocorrelation) and its
    best fitting exponential function

    :param ifg_path: An interferogram.
        ifg: :py:class:`pyrate.shared.Ifg`.
    :param: params: dict
        dict of config params
    :param calc_alpha: bool
        whether you calculate alpha.
    """
    # pylint: disable=invalid-name
    # pylint: disable=too-many-locals
    if isinstance(ifg_path, str):  # used during MPI
        ifg = shared.Ifg(ifg_path)
        ifg.open()
    else:
        ifg = ifg_path
    # assert isinstance(ifg_path, shared.Ifg)
    # ifg = ifg_path
    shared.nan_and_mm_convert(ifg, params)
    # calculate 2D auto-correlation of image using the
    # spectral method (Wiener-Khinchin theorem)
    if ifg.nan_converted:  # saves heaps of time with no-nan conversion
        phase = where(isnan(ifg.phase_data), 0, ifg.phase_data)
    else:
        phase = ifg.phase_data
    # distance division factor of 1000 converts to km and is needed to match
    # Matlab code output
    distfact = 1000

    nrows, ncols = phase.shape
    fft_phase = fft2(phase)
    pspec = real(fft_phase)**2 + imag(fft_phase)**2
    autocorr_grid = ifft2(pspec)
    nzc = np.sum(np.sum(phase != 0))
    autocorr_grid = fftshift(real(autocorr_grid)) / nzc

    # pixel distances from pixel at zero lag (image centre).
    xx, yy = meshgrid(range(ncols), range(nrows))

    # r_dist is distance from the center
    # doing np.divide and np.sqrt will improve performance as it keeps
    # calculations in the numpy land
    r_dist = np.divide(np.sqrt(((xx-ifg.x_centre) * ifg.x_size)**2 +
                               ((yy-ifg.y_centre) * ifg.y_size)**2), distfact)

    r_dist = reshape(r_dist, ifg.num_cells)
    acg = reshape(autocorr_grid, ifg.num_cells)

    # Symmetry in image; keep only unique points
    # tmp = unique_points(zip(acg, r_dist))
    # Sudipta: Is this faster than keeping only the 1st half as in Matlab?
    # Sudipta: Unlikely, as unique_point is a search/comparison,
    # whereas keeping 1st half is just numpy indexing.
    # If it is not faster, why was this done differently here?

    r_dist = r_dist[:int(ceil(ifg.num_cells/2.0)) + ifg.nrows]
    acg = acg[:len(r_dist)]

    # Alternative method to remove duplicate cells from Matlab Pirate
    # r_dist = r_dist[:ceil(len(r_dist)/2)+nlines]
    #  Reason for '+nlines' term unknown

    # eg. array([x for x in set([(1,1), (2,2), (1,1)])])
    # the above shortens r_dist by some number of cells

    # bin width for collecting data
    bin_width = max(ifg.x_size, ifg.y_size) * 2 / distfact

    # pick the smallest axis to determine circle search radius
    # print 'ifg.X_CENTRE, ifg.Y_CENTRE=', ifg.x_centre, ifg.y_centre
    # print 'ifg.X_SIZE, ifg.Y_SIZE', ifg.x_size, ifg.y_size
    if (ifg.x_centre * ifg.x_size) < (ifg.y_centre * ifg.y_size):
        maxdist = ifg.x_centre * ifg.x_size / distfact
    else:
        maxdist = ifg.y_centre * ifg.y_size/ distfact

    # filter out data where the of lag distance is greater than maxdist
    # r_dist = array([e for e in rorig if e <= maxdist]) #
    # MG: prefers to use all the data
    # acg = array([e for e in rorig if e <= maxdist])
    indices_to_keep = r_dist < maxdist
    r_dist = r_dist[indices_to_keep]
    acg = acg[indices_to_keep]

    if isinstance(ifg_path, str):
        ifg.close()

    if calc_alpha:
        # classify values of r_dist according to bin number
        rbin = ceil(r_dist / bin_width).astype(int)
        maxbin = max(rbin)  # consistent with Matlab code

        cvdav = zeros(shape=(2, maxbin))

        # the following stays in numpy land
        # distance instead of bin number
        cvdav[0, :] = np.multiply(range(maxbin), bin_width)
        # mean variance for the bins
        cvdav[1, :] = [mean(acg[rbin == b]) for b in range(maxbin)]

        # calculate best fit function maxvar*exp(-alpha*r_dist)
        alphaguess = 2 / (maxbin * bin_width)
        alpha = fmin(pendiffexp, x0=alphaguess, args=(cvdav,), disp=0,
                     xtol=1e-6, ftol=1e-6)
        print("1st guess alpha", alphaguess, 'converged alpha:', alpha)
        # maximum variance usually at the zero lag: max(acg[:len(r_dist)])
        return np.max(acg), alpha[0]
    else:
        return np.max(acg), None


def get_vcmt(ifgs, maxvar):
    """
    Returns the temporal variance/covariance matrix.
    """
    # pylint: disable=too-many-locals
    # c=0.5 for common master or slave; c=-0.5 if master
    # of one matches slave of another

    if isinstance(ifgs, dict):
        from collections import OrderedDict
        ifgs = {k: v for k, v in ifgs.items() if isinstance(v, PrereadIfg)}
        ifgs = OrderedDict(sorted(ifgs.items()))
        # pylint: disable=redefined-variable-type
        ifgs = ifgs.values()

    nifgs = len(ifgs)
    vcm_pat = zeros((nifgs, nifgs))

    dates = [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]
    ids = master_slave_ids(dates)

    for i, ifg in enumerate(ifgs):
        mas1, slv1 = ids[ifg.master], ids[ifg.slave]

        for j, ifg2 in enumerate(ifgs):
            mas2, slv2 = ids[ifg2.master], ids[ifg2.slave]
            if mas1 == mas2 or slv1 == slv2:
                vcm_pat[i, j] = 0.5

            if mas1 == slv2 or slv1 == mas2:
                vcm_pat[i, j] = -0.5

            if mas1 == mas2 and slv1 == slv2:
                vcm_pat[i, j] = 1.0  # handle testing ifg against itself

    # make covariance matrix in time domain
    std = sqrt(maxvar).reshape((nifgs, 1))
    vcm_t = std * std.transpose()
    return vcm_t * vcm_pat
