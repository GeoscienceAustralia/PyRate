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
This Python module implements a spatial low pass filter.
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import logging
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift
from scipy.interpolate import griddata
from pyrate import config as cf
from pyrate.covariance import cvd_from_phase, RDist

log = logging.getLogger(__name__)


def spatial_low_pass_filter(ts_lp, ifg, params):
    """
    Filter time series data spatially using either a Butterworth or Gaussian
    low pass filter defined by a cut-off distance. If the cut-off distance is
    defined as zero in the parameters dictionary then it is calculated for
    each time step using the pyrate.covariance.cvd_from_phase method.

    :param ndarray ts_lp: Array of time series data, the result of a temporal
                low pass filter operation. shape (ifg.shape, n_epochs)
    :param shared.Ifg instance ifg: interferogram object
    :param dict params: Dictionary of configuration parameters

    :return: ts_hp: filtered time series data of shape (ifg.shape, n_epochs)
    :rtype: ndarray
    """
    log.info('Applying spatial low pass filter')
    if params[cf.SLPF_NANFILL] == 0:
        ts_lp[np.isnan(ts_lp)] = 0  # need it here for cvd and fft
    else:  # optionally interpolate, operation is inplace
        _interpolate_nans(ts_lp, params[cf.SLPF_NANFILL_METHOD])
    r_dist = RDist(ifg)()
    for i in range(ts_lp.shape[2]):
        ts_lp[:, :, i] = _slpfilter(ts_lp[:, :, i], ifg, r_dist, params)
    log.info('Finished applying spatial low pass filter')
    return ts_lp


def _interpolate_nans(arr, method='linear'):
    """
    Fill any NaN values in arr with interpolated values. Nanfill and
    interpolation are performed in place.
    """
    rows, cols = np.indices(arr.shape[:2])
    for i in range(arr.shape[2]):
        a = arr[:, :, i]
        _interpolate_nans_2d(a, rows, cols, method)


def _interpolate_nans_2d(a, rows, cols, method):
    """
    In-place array interpolation and nanfill

    :param ndarray a: 2d ndarray to be interpolated
    :param ndarray rows: 2d ndarray of row indices
    :param ndarray cols: 2d ndarray of col indices
    :param str method: Method; one of 'nearest', 'linear', and 'cubic'
    """
    a[np.isnan(a)] = griddata(
        (rows[~np.isnan(a)], cols[~np.isnan(a)]),  # points we know
        a[~np.isnan(a)],  # values we know
        (rows[np.isnan(a)], cols[np.isnan(a)]),  # points to interpolate
        method=method
    )
    a[np.isnan(a)] = 0  # zero fill boundary/edge nans


def _slpfilter(phase, ifg, r_dist, params):
    """
    Wrapper function for spatial low pass filter
    """
    if np.all(np.isnan(phase)):  # return for nan matrix
        return phase
    cutoff = params[cf.SLPF_CUTOFF]

    if cutoff == 0:
        _, alpha = cvd_from_phase(phase, ifg, r_dist, calc_alpha=True)
        cutoff = 1.0/alpha
    rows, cols = ifg.shape
    return _slp_filter(phase, cutoff, rows, cols,
                       ifg.x_size, ifg.y_size, params)


def _slp_filter(phase, cutoff, rows, cols, x_size, y_size, params):
    """
    Function to perform spatial low pass filter
    """
    cx = np.floor(cols/2)
    cy = np.floor(rows/2)
    # fft for the input image
    imf = fftshift(fft2(phase))
    # calculate distance
    distfact = 1.0e3  # to convert into meters
    [xx, yy] = np.meshgrid(range(cols), range(rows))
    xx = (xx - cx) * x_size  # these are in meters as x_size in meters
    yy = (yy - cy) * y_size
    dist = np.sqrt(xx ** 2 + yy ** 2)/distfact  # km

    if params[cf.SLPF_METHOD] == 1:  # butterworth low pass filter
        H = 1. / (1 + ((dist / cutoff) ** (2 * params[cf.SLPF_ORDER])))
    else:  # Gaussian low pass filter
        H = np.exp(-(dist ** 2) / (2 * cutoff ** 2))
    outf = imf * H
    out = np.real(ifft2(ifftshift(outf)))
    out[np.isnan(phase)] = np.nan
    return out  # out is units of phase, i.e. mm
