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
This Python module implements pixel-by-pixel rate
(velocity) estimation using an iterative weighted least-squares
stacking method.
"""
import itertools

from scipy.linalg import solve, cholesky, qr, inv
from numpy import nan, isnan, sqrt, diag, delete, array, float32, size
import numpy as np
from joblib import Parallel, delayed
from pyrate.core import config as cf
from pyrate.core.shared import joblib_log_level
from pyrate.core.logger import pyratelogger as log


def stack_rate_array(ifgs, params, vcmt, mst=None):
    """
    This function loops over all interferogram pixels in a 3-dimensional array and estimates
    the pixel rate (velocity) by applying the iterative weighted least-squares stacking
    algorithm 'pyrate.core.stack.stack_rate_pixel'.

    :param Ifg.object ifgs: Sequence of objects containing the interferometric observations
    :param dict params: Configuration parameters
    :param ndarray vcmt: Derived positive definite temporal variance covariance matrix
    :param ndarray mst: Pixel-wise matrix describing the minimum spanning tree network

    :return: rate: Rate (velocity) map
    :rtype: ndarray
    :return: error: Standard deviation of the rate map
    :rtype: ndarray
    :return: samples: Number of observations used in rate calculation for each pixel
    :rtype: ndarray
    """
    nsig, pthresh, cols, error, mst, obs, parallel, _, rate, rows, samples, span = _stack_setup(ifgs, mst, params)

    # pixel-by-pixel calculation.
    # nested loops to loop over the 2 image dimensions
    if parallel:
        log.info('Calculating stack rate in parallel')
        res = Parallel(n_jobs=params[cf.PROCESSES], verbose=joblib_log_level(cf.LOG_LEVEL))(
            delayed(stack_rate_pixel)(obs[:, r, c], mst[:, r, c], vcmt, span, nsig, pthresh) for r, c in itertools.product(range(rows), range(cols))
        )
        res = np.array(res)

        rate = res[:, 0].reshape(rows, cols)
        error = res[:, 1].reshape(rows, cols)
        samples = res[:, 2].reshape(rows, cols)
    else:
        log.info('Calculating stack rate in serial')
        for i in range(rows):
            for j in range(cols):
                rate[i, j], error[i, j], samples[i, j] = stack_rate_pixel(obs[:, i, j], mst[:, i, j], vcmt, span, nsig, pthresh)

    return rate, error, samples


def mask_rate(rate, error, maxsig):
    """
    Function to mask pixels in the rate and error arrays when the error
    is greater than the error threshold 'maxsig'.

    :param ndarray rate: array of pixel rates derived by stacking
    :param ndarray error: array of errors for the pixel rates
    :param int maxsig: error threshold for masking (in millimetres).

    :return: rate: Masked rate (velocity) map
    :rtype: ndarray
    :return: error: Masked error (standard deviation) map
    :rtype: ndarray
    """
    log.info('Masking stack rate and error maps where sigma is greater than {} millimetres'.format(maxsig))
    # initialise mask array with existing NaNs
    mask = ~isnan(error)
    # original Nan count
    orig = np.count_nonzero(mask)
    # determine where error is larger than the maximum sigma threshold
    mask[mask] &= error[mask] > maxsig
    # replace values with NaNs
    rate[mask] = nan
    error[mask] = nan
    # calculate percentage of masked pixels
    nummasked = int(np.count_nonzero(mask)/orig*100)
    log.info('Percentage of pixels masked = {}%'.format(nummasked))

    return rate, error


def stack_rate_pixel(obs, mst, vcmt, span, nsig, pthresh):
    """
    Algorithm to estimate the rate (velocity) for a single pixel using iterative
    weighted least-squares stacking method.

    :param ndarray obs: Vector of interferometric phase observations for the pixel
    :param ndarray mst: Vector describing the minimum spanning tree network for the pixel
    :param ndarray vcmt: Derived positive definite temporal variance covariance matrix
    :param ndarray span: Vector of interferometric time spans
    :param int nsig: Threshold for iterative removal of interferometric observations
    :param int pthresh: Threshold for minimum number of observations for the pixel

    :return: rate: Estimated rate (velocity) for the pixel
    :rtype: float64
    :return: error: Standard deviation of the observations with respect to the estimated rate
    :rtype: float64
    :return: samples: Number of observations used in the rate estimation for the pixel
    :rtype: int
    """

    # find the indices of independent ifgs from MST
    ind = np.nonzero(mst)[0]  # only True's in mst are chosen
    # iterative loop to calculate 'robust' velocity for pixel
    default_no_samples = len(ind)

    while len(ind) >= pthresh:
        # select ifg observations
        ifgv = obs[ind]

        # form design matrix from appropriate ifg time spans
        B = span[:, ind]

        # Subset of full VCM matrix for selected observations
        vcm_temp = vcmt[ind, np.vstack(ind)]

        # Get the lower triangle cholesky decomposition.
        # V must be positive definite (symmetrical and square)
        T = cholesky(vcm_temp, 1)

        # Incorporate inverse of VCM into the design matrix
        # and observations vector
        A = solve(T, B.transpose())
        b = solve(T, ifgv.transpose())

        # Factor the design matrix, incorporate covariances or weights into the
        # system of equations, and transform the response vector.
        Q, R, _ = qr(A, mode='economic', pivoting=True)
        z = Q.conj().transpose().dot(b)

        # Compute the Lstsq coefficient for the velocity
        v = solve(R, z)

        # Compute the model errors
        err1 = inv(vcm_temp).dot(B.conj().transpose())
        err2 = B.dot(err1)
        err = sqrt(diag(inv(err2)))

        # Compute the residuals (model minus observations)
        r = (B * v) - ifgv

        # determine the ratio of residuals and apriori variances
        w = cholesky(inv(vcm_temp))
        wr = abs(np.dot(w, r.transpose()))

        # test if maximum ratio is greater than user threshold.
        max_val = wr.max()
        if max_val > nsig:
            # if yes, discard and re-do the calculation.
            ind = delete(ind, wr.argmax())
        else:
            # if no, save estimate, exit the while loop and go to next pixel
            return v[0], err[0], ifgv.shape[0]
    # dummy return for no change
    return np.nan, np.nan, default_no_samples


def _stack_setup(ifgs, mst, params):
    """
    Convenience function for stack rate setup
    """
    # MULTIPROCESSING parameters
    parallel = params[cf.PARALLEL]
    processes = params[cf.PROCESSES]
    # stack rate parameters from config file
    # n-sigma ratio used to threshold 'model minus observation' residuals
    nsig = params[cf.LR_NSIG]
    # Pixel threshold; minimum number of observations for a pixel
    pthresh = params[cf.LR_PTHRESH]
    rows, cols = ifgs[0].phase_data.shape
    # make 3D block of observations
    obs = array([np.where(isnan(x.phase_data), 0, x.phase_data) for x in ifgs])
    span = array([[x.time_span for x in ifgs]])
    # Update MST in case additional NaNs generated by APS filtering
    if mst is None:  # dummy mst if none is passed in
        mst = ~isnan(obs)
    else:
        mst[isnan(obs)] = 0

    # preallocate empty arrays (no need to preallocate NaNs)
    error = np.empty([rows, cols], dtype=float32)
    rate = np.empty([rows, cols], dtype=float32)
    samples = np.empty([rows, cols], dtype=np.float32)
    return nsig, pthresh, cols, error, mst, obs, parallel, processes, rate, rows, samples, span
