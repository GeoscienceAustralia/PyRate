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
This Python module contains functions for computing a time series
inversion in PyRate. It is based on the Matlab codes 'tsinvnosm.m'
and 'tsinvlap.m' contained in the Pirate package.
"""
# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
import itertools
from numpy import (where, isnan, nan, diff, zeros,
                   float32, cumsum, dot, delete, asarray)
from numpy.linalg import matrix_rank, pinv, cholesky
import numpy as np
from scipy.linalg import qr
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

from pyrate.algorithm import master_slave_ids, get_epochs
from pyrate import config as cf
from pyrate.config import ConfigException
from pyrate import mst as mst_module


def time_series_setup(ifgs, mst, params):
    """ Convenience function setting up time series computation parameters"""
    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    # if mst is not a single tree then do interpolation
    if mst_module.mst_from_ifgs(ifgs)[1]:
        interp = 0
    else:
        interp = 1

    # Parallel Processing parameters
    parallel = params[cf.PARALLEL]
    # Time Series parameters
    tsmethod = params[cf.TIME_SERIES_METHOD]
    
    if tsmethod == 1 and params[cf.TIME_SERIES_SM_ORDER] is None:
        missing_option_error(cf.TIME_SERIES_SM_ORDER)
    else:
        smorder = params[cf.TIME_SERIES_SM_ORDER]

    if tsmethod == 1 and params[cf.TIME_SERIES_SM_FACTOR] is None:
        missing_option_error(cf.TIME_SERIES_SM_FACTOR)
    else:
        smfactor = np.power(10, params[cf.TIME_SERIES_SM_FACTOR])   

    if params[cf.TIME_SERIES_PTHRESH] is None:
        missing_option_error(cf.TIME_SERIES_PTHRESH)
    else:
        pthresh = params[cf.TIME_SERIES_PTHRESH]
        
    if pthresh < 0.0 or pthresh > 1000:
        raise ValueError(
            "minimum number of coherent observations for a pixel"
            " TIME_SERIES_PTHRESH setting must be >= 0.0 and <= 1000")

    epochlist = get_epochs(ifgs)[0]
    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    span = diff(epochlist.spans)
    nepoch = len(epochlist.dates)  # epoch number
    nvelpar = nepoch - 1  # velocity parameters number
    # nlap = nvelpar - smorder  # Laplacian observations number
    mast_slave_ids = master_slave_ids(epochlist.dates)
    imaster = [mast_slave_ids[ifg.master] for ifg in ifgs]
    islave = [mast_slave_ids[ifg.slave] for ifg in ifgs]
    imaster = min(imaster, islave)
    islave = max(imaster, islave)
    b0_mat = zeros((nifgs, nvelpar))
    for i in range(nifgs):
        b0_mat[i, imaster[i]:islave[i]] = span[imaster[i]:islave[i]]

    # change the sign if slave is earlier than master
    isign = where(imaster > islave)
    b0_mat[isign[0], :] = -b0_mat[isign[0], :]
    tsvel_matrix = np.empty(shape=(nrows, ncols, nvelpar),
                            dtype=float32)
    ifg_data = np.zeros((nifgs, nrows, ncols), dtype=float32)
    for ifg_num in range(nifgs):
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data
    if mst is None:
        mst = ~isnan(ifg_data)
    return b0_mat, interp, pthresh, smfactor, smorder, tsmethod, ifg_data, \
        mst, ncols, nrows, nvelpar, parallel, span, tsvel_matrix


def time_series(ifgs, params, vcmt, mst=None):
    """
    Returns time series data from the given interferograms.

    :param ifgs: Network of interferograms.
    :param params: configuration parameters
    :param vcmt: Derived positive definite temporal variance covariance matrix
    :param mst: [optional] array of ifg indexes from the MST-matrix.
    :param parallel: use parallel processing or not.

    :return: Tuple with the elements:

        - incremental displacement time series,
        - cumulative displacement time series, and
        - velocity for the epoch interval

        these outputs are multi-dimensional arrays of
        size(nrows, ncols, nepochs-1), where:

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and
        - *nepochs* is the number of unique epochs (dates)
        covered by the ifgs.).
    """

    b0_mat, interp, p_thresh, sm_factor, sm_order, ts_method, ifg_data, mst, \
        ncols, nrows, nvelpar, parallel, span, tsvel_matrix = \
        time_series_setup(ifgs, mst, params)

    if parallel == 1:
        tsvel_matrix = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(time_series_by_rows)(r, b0_mat, sm_factor, sm_order,
                                         ifg_data, mst, ncols, nvelpar,
                                         p_thresh, vcmt, ts_method, interp)
            for r in range(nrows))

    elif parallel == 2:

        res = np.array(Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(time_series_by_pixel)(i, j, b0_mat, sm_factor, sm_order,
                                          ifg_data, mst, nvelpar, p_thresh,
                                          vcmt, ts_method, interp)
            for (i, j) in itertools.product(range(nrows), range(ncols))))
        tsvel_matrix = np.reshape(res, newshape=(nrows, ncols, res.shape[1]))
    else:
        for row in range(nrows):
            for col in range(ncols):
                tsvel_matrix[row, col] = time_series_by_pixel(
                    row, col, b0_mat, sm_factor, sm_order, ifg_data, mst,
                    nvelpar, p_thresh, vcmt, ts_method, interp)

    tsvel_matrix = where(tsvel_matrix == 0, nan, tsvel_matrix)
    # SB: do the span multiplication as a numpy linalg operation, MUCH faster
    #  not even this is necessary here, perform late for performance
    tsincr = tsvel_matrix * span
    tscum = cumsum(tsincr, 2)
    # SB: convert zeros to nans not performed in matlab
    # SB: perform this after tsvel_matrix has been nan converted,
    # saves the step of comparing a large matrix (tsincr) to zero.
    # tscum = where(tscum == 0, nan, tscum)
    return tsincr, tscum, tsvel_matrix


def time_series_by_rows(row, b0_mat, sm_factor, sm_order, ifg_data, mst, ncols,
                        nvelpar, p_thresh, vcmt, ts_method, interp):
    """ time series computation for each row of interferrograms """
    tsvel = np.empty(shape=(ncols, nvelpar), dtype=float32)
    for col in range(ncols):
        tsvel[col, :] = time_series_by_pixel(
            row, col, b0_mat, sm_factor, sm_order, ifg_data, mst, nvelpar,
            p_thresh, vcmt, ts_method, interp)

    return tsvel


def remove_rank_def_rows(b_mat, nvelpar, ifgv, sel):
    """ Remove rank deficient rows """
    _, _, e_var = qr(b_mat, mode='economic', pivoting=True)
    licols = e_var[matrix_rank(b_mat):nvelpar]
    [rmrow, _] = where(b_mat[:, licols] != 0)
    b_mat = delete(b_mat, rmrow, axis=0)
    ifgv = delete(ifgv, rmrow)
    sel = delete(sel, rmrow)
    return b_mat, ifgv, sel, rmrow


def time_series_by_pixel(row, col, b0_mat, sm_factor, sm_order, ifg_data, mst,
                         nvelpar, p_thresh, vcmt, method, interp):
    """ time series computation for each pixel """
    # check pixel for non-redundant ifgs
    sel = np.nonzero(mst[:, row, col])[0]  # trues in mst are chosen
    if len(sel) >= p_thresh:
        ifgv = ifg_data[sel, row, col]
        # make design matrix, b_mat
        b_mat = b0_mat[sel, :]
        if interp == 0:
            # remove rank deficient rows
            rmrow = asarray([0])  # dummy

            while len(rmrow) > 0:
                # if b_mat.shape[0] <=1 then we return nans
                if b_mat.shape[0] > 1:
                    b_mat, ifgv, sel, rmrow = remove_rank_def_rows(
                        b_mat, nvelpar, ifgv, sel)
                else:
                    return np.empty(nvelpar) * np.nan

            # Some epochs have been deleted; get valid epoch indices
            velflag = sum(abs(b_mat), 0)
            # remove corresponding columns in design matrix
            b_mat = b_mat[:, ~np.isclose(velflag, 0.0)]
        else:
            velflag = np.ones(nvelpar)
        if method == 1:
            # Use Laplacian smoothing method
            tsvel = _solve_ts_lap(nvelpar, velflag, ifgv, b_mat,
                                  sm_order, sm_factor, sel, vcmt)
        elif method == 2:
            # Use SVD method
            tsvel = solve_ts_svd(nvelpar, velflag, ifgv, b_mat)
        else:
            raise ValueError("Unrecognised time series method")
        return tsvel
    else:
        return np.empty(nvelpar) * np.nan


def solve_ts_svd(nvelpar, velflag, ifgv, b_mat):
    """
    Solve the linear least squares system using the SVD method.
    Similar to the SBAS method implemented by Berardino et al. 2002
    """
    # pre-allocate the velocity matrix
    tsvel = np.empty(nvelpar, dtype=float32) * np.nan
    # solve least squares equation using Moore-Penrose pseudoinverse
    tsvel[velflag != 0] = dot(pinv(b_mat), ifgv)
    return tsvel


def _solve_ts_lap(nvelpar, velflag, ifgv, mat_b, smorder, smfactor, sel, vcmt):
    """
    Solve the linear least squares system using the Finite Difference
    method using a Laplacian Smoothing operator.
    Similar to the method implemented by Schmidt and Burgmann 2003
    """
    # pylint: disable=invalid-name
    # Laplacian observations number
    nlap = nvelpar - smorder
    #  Laplacian smoothing coefficient matrix
    b_lap0 = np.zeros(shape=(nlap, nvelpar))

    for i in range(nlap):
        if smorder == 1:
            b_lap0[i, i:i+2] = [-1, 1]
        else:
            b_lap0[i, i:i+3] = [1, -2, 1]

    # Scale the coefficients by Laplacian smoothing factor
    b_lap0 *= smfactor

    # Laplacian smoothing design matrix
    nvelleft = np.count_nonzero(velflag)
    nlap = nvelleft - smorder

    # constrain for the first and the last incremental
    b_lap1 = - np.divide(np.ones(shape=nvelleft), nvelleft - 1)
    b_lap1[0] = 1.0
    b_lapn = - np.divide(np.ones(shape=nvelleft), nvelleft - 1)
    b_lapn[-1] = 1.0

    b_lap = np.empty(shape=(nlap + 2, nvelleft))
    b_lap[0, :] = b_lap1
    b_lap[1:nlap + 1, :] = b_lap0[0:nlap, 0:nvelleft]
    b_lap[-1, :] = b_lapn

    nlap += 2

    # add laplacian design matrix to existing design matrix
    mat_b = np.concatenate((mat_b, b_lap), axis=0)

    # combine ifg and Laplacian smooth vector
    v_lap = np.zeros(nlap)
    obsv = np.concatenate((ifgv, v_lap), axis=0)

    # make variance-covariance matrix
    # new covariance matrix, adding the laplacian equations
    m = len(sel)
    nobs = m + nlap
    vcm_tmp = np.eye(nobs)
    vcm_tmp[:m, :m] = vcmt[sel, np.vstack(sel)]

    # solve the equation by least-squares
    # calculate velocities
    # we get the lower triangle in numpy, matlab gives upper triangle
    w = cholesky(pinv(vcm_tmp)).T
    wb = dot(w, mat_b)
    wl = dot(w, obsv)
    x = dot(pinv(wb, rcond=1e-8), wl)

    # TODO: implement residuals and roughness calculations
    tsvel = np.empty(nvelpar, dtype=float32) * np.nan
    tsvel[~np.isclose(velflag, 0.0, atol=1e-8)] = x[:nvelleft]

    # TODO: implement uncertainty estimates (tserror) like in matlab code
    return tsvel


def plot_timeseries(tsincr, tscum, tsvel, output_dir):  # pragma: no cover
    """
    Plot the results of the timeseries analysis.

    A very basic plotting function used in code development

    .. todo:: this should be moved and tidied up
    """

    nvelpar = len(tsincr[0, 0, :])
    for i in range(nvelpar):

        imgplot = plt.imshow(tsincr[:, :, i])
        imgplot.set_clim(-10, 10)
        plt.colorbar(ticks=[-10, 0, 10], orientation='horizontal')
        plt.draw()
        plt.savefig(output_dir + 'tsincr_' + str(i) + '.png')
        plt.close()
        plt.figure()
        imgplot = plt.imshow(tscum[:, :, i])
        imgplot.set_clim(-10, 10)
        plt.colorbar(ticks=[-10, 0, 10], orientation='horizontal')
        plt.draw()
        plt.savefig(output_dir + 'tscum_' + str(i) + '.png')
        plt.close()

        imgplot = plt.imshow(tsvel[:, :, i])
        imgplot.set_clim(-50, 100)
        plt.colorbar(ticks=[-50, 50, 100], orientation='horizontal')
        plt.draw()
        plt.savefig(output_dir + 'tsvel_' + str(i) + '.png')
        plt.close()


def missing_option_error(option):
    '''
    Internal convenience function for raising similar missing option errors.
    '''
    msg = "Missing '%s' option in config file" % option
    raise ConfigException(msg)

class TimeSeriesError(Exception):
    """
    Generic exception for time series errors.
    """
