#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
Variance/Covariance matrix functionality.
"""
# coding: utf-8
from os.path import basename, join
from collections import OrderedDict
from numpy import array, where, isnan, real, imag, sqrt, meshgrid
from numpy import zeros, vstack, ceil, mean, exp, reshape
from numpy.linalg import norm
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.optimize import fmin

import pyrate.constants as C
from pyrate.core import shared, ifgconstants as ifc, mpiops
from pyrate.core.shared import PrereadIfg, Ifg
from pyrate.core.algorithm import first_second_ids
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import Configuration

# pylint: disable=too-many-arguments
MAIN_PROCESS = 0


def _pendiffexp(alphamod, cvdav):
    """
    Fits an exponential model to data.

    :param float alphamod: Exponential decay exponent.
    :param ndarray cvdav: Function magnitude at 0 radius (2 col array of
    radius, variance)
    """
    # pylint: disable=invalid-name
    # maxvar usually at zero lag
    mx = cvdav[1, 0]
    return norm(cvdav[1, :] - (mx * exp(-alphamod * cvdav[0, :])))


# this is not used any more
def _unique_points(points):  # pragma: no cover
    """
    Returns unique points from a list of coordinates.

    :param points: Sequence of (y,x) or (x,y) tuples.
    """
    return vstack([array(u) for u in set(points)])


def cvd(ifg_path, params, r_dist, calc_alpha=False, write_vals=False, save_acg=False):
    """
    Calculate the 1D covariance function of an entire interferogram as the
    radial average of its 2D autocorrelation.

    :param str ifg_path: An interferogram file path. OR
    :param dict params: Dictionary of configuration parameters
    :param ndarray r_dist: Array of distance values from the image centre
                (See Rdist class for more details)
    :param bool calc_alpha: If True calculate alpha
    :param bool write_vals: If True write maxvar and alpha values to
                interferogram metadata
    :param bool save_acg: If True write autocorrelation and radial distance
                data to numpy array file on disk

    :return: maxvar: The maximum variance (at zero lag)
    :rtype: float
    :return: alpha: the exponential length-scale of decay factor
    :rtype: float
    """
    ifg = shared.Ifg(ifg_path)
    ifg.open()
    shared.nan_and_mm_convert(ifg, params)
    # calculate 2D auto-correlation of image using the
    # spectral method (Wiener-Khinchin theorem)
    if ifg.nan_converted:  # if nancoverted earlier, convert nans back to 0's
        phase = where(isnan(ifg.phase_data), 0, ifg.phase_data)
    else:
        phase = ifg.phase_data

    maxvar, alpha = cvd_from_phase(phase, ifg, r_dist, calc_alpha, save_acg=save_acg, params=params)
    if write_vals:
        ifg.add_metadata(**{
            ifc.PYRATE_MAXVAR: str(maxvar),
            ifc.PYRATE_ALPHA: str(alpha)
        })

    ifg.close()

    return maxvar, alpha


def _save_cvd_data(acg, r_dist, ifg_path, outdir):
    """
    Function to save numpy array of autocorrelation data to disk
    """
    data = np.column_stack((acg, r_dist))
    suffix = basename(ifg_path).split('.')[0]
    data_file = join(outdir, f'cvd_data_{suffix}.npy')
    np.save(file=data_file, arr=data)


def cvd_from_phase(phase, ifg, r_dist, calc_alpha, save_acg=False, params=None):
    """
    A convenience function used to compute radial autocovariance from phase
    data

    :param ndarray phase: An array of interferogram phase data
    :param Ifg class ifg: A pyrate.shared.Ifg class instance
    :param ndarray r_dist: Array of distance values from the image centre
                (See Rdist class for more details)
    :param bool calc_alpha: If True calculate alpha
    :param bool save_acg: If True write autocorrelation and radial distance
                data to numpy array file on disk
    :param dict params: [optional] Dictionary of configuration parameters;
                Must be provided if save_acg=True

    :return: maxvar: The maximum variance (at zero lag)
    :rtype: float
    :return: alpha: the exponential length-scale of decay factor
    :rtype: float
    """
    # pylint: disable=invalid-name
    # pylint: disable=too-many-locals

    autocorr_grid = _get_autogrid(phase)
    acg = reshape(autocorr_grid, phase.size, order='F')
    # Symmetry in image; keep only unique points
    # tmp = _unique_points(zip(acg, r_dist))
    # Sudipta: Unlikely, as unique_point is a search/comparison,
    # whereas keeping 1st half is just numpy indexing.
    # If it is not faster, why was this done differently here?
    # r_dist = r_dist[:int(ceil(phase.size / 2.0)) + nrows]
    acg = acg[:len(r_dist)]
    # Alternative method to remove duplicate cells
    # r_dist = r_dist[:ceil(len(r_dist)/2)+nlines]
    #  Reason for '+nlines' term unknown
    # eg. array([x for x in set([(1,1), (2,2), (1,1)])])
    # the above shortens r_dist by some number of cells

    # pick the smallest axis to determine circle search radius
    if (ifg.x_centre * ifg.x_size) < (ifg.y_centre * ifg.y_size):
        maxdist = (ifg.x_centre+1) * ifg.x_size / ifc.METRE_PER_KM
    else:
        maxdist = (ifg.y_centre+1) * ifg.y_size / ifc.METRE_PER_KM

    # filter out data where the of lag distance is greater than maxdist
    # r_dist = array([e for e in rorig if e <= maxdist]) #
    # MG: prefers to use all the data
    # acg = array([e for e in rorig if e <= maxdist])
    indices_to_keep = r_dist < maxdist
    acg = acg[indices_to_keep]

    # optionally save acg vs dist observations to disk
    if save_acg:
        _save_cvd_data(acg, r_dist[indices_to_keep],
                       ifg.data_path, params[C.TMPDIR])

    if calc_alpha:
        # bin width for collecting data
        bin_width = max(ifg.x_size, ifg.y_size) * 2 / ifc.METRE_PER_KM  # km
        r_dist = r_dist[indices_to_keep]  # km
        # classify values of r_dist according to bin number
        rbin = ceil(r_dist / bin_width).astype(int)
        maxbin = max(rbin) - 1  # consistent with Legacy data

        cvdav = zeros(shape=(2, maxbin + 1))

        # the following stays in numpy land
        # distance instead of bin number
        cvdav[0, :] = np.multiply(range(maxbin + 1), bin_width)
        # mean variance for the bins
        cvdav[1, :] = [mean(acg[rbin == b]) for b in range(maxbin + 1)]
        # calculate best fit function maxvar*exp(-alpha*r_dist)
        alphaguess = 2 / (maxbin * bin_width)
        alpha = fmin(_pendiffexp, x0=alphaguess, args=(cvdav,), disp=False,
                     xtol=1e-6, ftol=1e-6)
        log.debug(f"1st guess alpha {alphaguess}, converged alpha: {alpha}")
        # maximum variance usually at the zero lag: max(acg[:len(r_dist)])
        return np.max(acg), alpha[0]  # alpha unit 1/km

    return np.max(acg), None


class RDist():
    """
    RDist class used for caching r_dist during maxvar/alpha computation
    """
    # pylint: disable=invalid-name
    def __init__(self, ifg):
        self.r_dist = None
        self.ifg = ifg
        self.nrows, self.ncols = ifg.shape

    def __call__(self):

        if self.r_dist is None:
            size = self.nrows * self.ncols
            # pixel distances from pixel at zero lag (image centre).
            xx, yy = meshgrid(range(self.ncols), range(self.nrows))
            # r_dist is distance from the center
            # doing np.divide and np.sqrt will improve performance as it keeps
            # calculations in the numpy land
            self.r_dist = np.divide(np.sqrt(((xx - self.ifg.x_centre) *
                                             self.ifg.x_size) ** 2 +
                                            ((yy - self.ifg.y_centre) *
                                             self.ifg.y_size) ** 2),
                                    ifc.METRE_PER_KM)  # km
            self.r_dist = reshape(self.r_dist, size, order='F')
            self.r_dist = self.r_dist[:int(ceil(size / 2.0)) + self.nrows]

        return self.r_dist


def _get_autogrid(phase):
    """
    Helper function to assist with memory re-allocation during FFT calculation
    """
    autocorr_grid = _calc_autoc_grid(phase)
    nzc = np.sum(np.sum(phase != 0))
    autocorr_grid = fftshift(real(autocorr_grid)) / nzc
    return autocorr_grid


def _calc_autoc_grid(phase):
    """
    Helper function to assist with memory re-allocation during FFT calculation
    """
    pspec = _calc_power_spectrum(phase)
    autocorr_grid = ifft2(pspec)
    return autocorr_grid.astype(dtype=np.complex64)


def _calc_power_spectrum(phase):
    """
    Helper function to assist with memory re-allocation during FFT calculation
    """
    fft_phase = fft2(phase)
    pspec = real(fft_phase) ** 2 + imag(fft_phase) ** 2
    return pspec.astype(dtype=np.float32)


def get_vcmt(ifgs, maxvar):
    """
    Assembles a temporal variance/covariance matrix using the method
    described by Biggs et al., Geophys. J. Int, 2007. Matrix elements are
    evaluated according to sig_i * sig_j * C_ij where i and j are two
    interferograms and C is a matrix of coefficients:

    C = 1 if the first and second epochs of i and j are equal
    C = 0.5 if have i and j share either a common first or second epoch
    C = -0.5 if the first of i or j equals the second of the other
    C = 0 otherwise

    :param list ifgs: A list of pyrate.shared.Ifg class objects.
    :param ndarray maxvar: numpy array of maximum variance values for the
                interferograms.

    :return: vcm_t: temporal variance-covariance matrix
    :rtype: ndarray
    """
    # pylint: disable=too-many-locals
    # c=0.5 for common first or second; c=-0.5 if first
    # of one matches second of another

    if isinstance(ifgs, dict):
        ifgs = {k: v for k, v in ifgs.items() if isinstance(v, PrereadIfg)}
        ifgs = OrderedDict(sorted(ifgs.items()))
        ifgs = ifgs.values()

    nifgs = len(ifgs)
    vcm_pat = zeros((nifgs, nifgs))

    dates = [ifg.first for ifg in ifgs] + [ifg.second for ifg in ifgs]
    ids = first_second_ids(dates)

    for i, ifg in enumerate(ifgs):
        mas1, slv1 = ids[ifg.first], ids[ifg.second]

        for j, ifg2 in enumerate(ifgs):
            mas2, slv2 = ids[ifg2.first], ids[ifg2.second]
            if mas1 == mas2 or slv1 == slv2:
                vcm_pat[i, j] = 0.5

            if mas1 == slv2 or slv1 == mas2:
                vcm_pat[i, j] = -0.5

            if mas1 == mas2 and slv1 == slv2:
                vcm_pat[i, j] = 1.0  # diagonal elements

    # make covariance matrix in time domain
    std = sqrt(maxvar).reshape((nifgs, 1))
    vcm_t = std * std.transpose()
    return vcm_t * vcm_pat


def maxvar_vcm_calc_wrapper(params):
    """
    MPI wrapper for maxvar and vcmt computation
    """
    preread_ifgs = params[C.PREREAD_IFGS]
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    log.info('Calculating the temporal variance-covariance matrix')

    def _get_r_dist(ifg_path):
        """
        Get RDIst class object
        """
        ifg = Ifg(ifg_path)
        ifg.open()
        r_dist = RDist(ifg)()
        ifg.close()
        return r_dist

    r_dist = mpiops.run_once(_get_r_dist, ifg_paths[0])
    prcs_ifgs = mpiops.array_split(list(enumerate(ifg_paths)))
    process_maxvar = {}
    n_prcs = len(prcs_ifgs)
    n_ifgs = len(ifg_paths)

    for n, i in prcs_ifgs:
        log.debug(f'Calculating maxvar for {n} of process ifgs {n_prcs} of total {n_ifgs}')
        val = cvd(i, params, r_dist, calc_alpha=True, write_vals=True, save_acg=True)[0]
        process_maxvar[int(n)] = val

    maxvar_d = shared.join_dicts(mpiops.comm.allgather(process_maxvar))
    maxvar = [v[1] for v in sorted(maxvar_d.items(), key=lambda s: s[0])]

    vcmt = mpiops.run_once(get_vcmt, preread_ifgs, maxvar)
    log.debug("Finished maxvar and vcm calc!")
    params[C.MAXVAR], params[C.VCMT] = maxvar, vcmt
    np.save(Configuration.vcmt_path(params), arr=vcmt)
    return maxvar, vcmt
