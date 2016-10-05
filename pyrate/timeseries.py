"""
Functions for doing an InSAR time series inversion in PyRate.
Based on the Matlab Pirate 'tsinvnosm.m' and 'tsinvlap.m' functions.

..codeauthor:: Vanessa Newey, Sudipta Basak, Matt Garthwaite
"""

from numpy import where, isnan, nan, diff, zeros, empty
from numpy import float32, cumsum, dot, delete, asarray
from numpy.linalg import matrix_rank, pinv, cholesky
from scipy.linalg import qr
import numpy as np
import matplotlib.pyplot as plt
import itertools
import parmap

import pyrate.config as config
from algorithm import master_slave_ids, get_epochs
import pyrate.ifgconstants as ifc
from pyrate import config as cf
from pyrate import mst as mst_module


def time_series_setup(ifgs, mst, params):
    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    # check if mst is not a tree, then do interpolate
    if mst_module.mst_from_ifgs(ifgs)[1]:
        params[cf.TIME_SERIES_INTERP] = 0
    else:
        params[cf.TIME_SERIES_INTERP] = 1

    # make sure by now we have the time series interpolation parameter
    assert params[cf.TIME_SERIES_INTERP] is not None

    # Parallel Processing parameters
    parallel = params[cf.PARALLEL]
    processes = params[cf.PROCESSES]
    # Time Series parameters
    TSMETHOD = params[cf.TIME_SERIES_METHOD]
    INTERP = params[cf.TIME_SERIES_INTERP]
    # INTERP = 1
    SMORDER = params[cf.TIME_SERIES_SM_ORDER]
    SMFACTOR = np.power(10, params[cf.TIME_SERIES_SM_FACTOR])
    PTHRESH = params[cf.TIME_SERIES_PTHRESH]
    head = ifgs[0]
    check_time_series_params(head, PTHRESH)
    epochlist = get_epochs(ifgs)
    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    span = diff(epochlist.spans)
    nepoch = len(epochlist.dates)  # epoch number
    nvelpar = nepoch - 1  # velocity parameters number
    nlap = nvelpar - SMORDER  # Laplacian observations number
    mast_slave_ids = master_slave_ids(epochlist.dates)
    imaster = [mast_slave_ids[ifg.master] for ifg in ifgs]
    islave = [mast_slave_ids[ifg.slave] for ifg in ifgs]
    imaster = min(imaster, islave)
    islave = max(imaster, islave)
    B0 = zeros((nifgs, nvelpar))
    for i in range(nifgs):
        B0[i, imaster[i]:islave[i]] = span[imaster[i]:islave[i]]

    # change the sign if slave is earlier than master
    isign = where(imaster > islave)
    B0[isign[0], :] = -B0[isign[0], :]
    tsvel_matrix = np.empty(shape=(nrows, ncols, nvelpar),
                            dtype=np.float32)
    ifg_data = np.zeros((nifgs, nrows, ncols), dtype=float32)
    for ifg_num in xrange(nifgs):
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data
    if mst is None:
        mst = ~isnan(ifg_data)
    return B0, INTERP, PTHRESH, SMFACTOR, SMORDER, TSMETHOD, ifg_data, mst, \
           ncols, nrows, nvelpar, parallel, processes, span, tsvel_matrix


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

        these outputs are multi-dimensional arrays of size(nrows, ncols, nepochs-1),
        where:

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and
        - *nepochs* is the number of unique epochs (dates) covered by the ifgs.).
    """

    B0, INTERP, PTHRESH, SMFACTOR, SMORDER, TSMETHOD, ifg_data, mst, \
    ncols, nrows, nvelpar, parallel, processes, span, tsvel_matrix = \
        time_series_setup(ifgs, mst, params)

    if parallel == 1:
        tsvel_matrix = parmap.map(time_series_by_rows, range(nrows), B0,
                                  SMFACTOR, SMORDER, ifg_data, mst, ncols,
                                  nvelpar, PTHRESH, vcmt, TSMETHOD, INTERP,
                                  processes=processes)
    elif parallel == 2:
        res = parmap.starmap(time_series_by_pixel,
                             itertools.product(range(nrows), range(ncols)),
                             B0, SMFACTOR, SMORDER, ifg_data, mst, nvelpar,
                             PTHRESH, vcmt, TSMETHOD, INTERP,
                             processes=processes)
        res = np.array(res)
        tsvel_matrix = np.reshape(res, newshape=(nrows, ncols, res.shape[1]))
    else:
        for row in xrange(nrows):
            for col in xrange(ncols):
                tsvel_matrix[row, col] = time_series_by_pixel(
                    row, col, B0, SMFACTOR, SMORDER, ifg_data, mst, nvelpar,
                    PTHRESH, vcmt, TSMETHOD, INTERP)

    tsvel_matrix = where(tsvel_matrix == 0, nan, tsvel_matrix)
    # SB: do the span multiplication as a numpy linalg operation, MUCH faster
    # not even this is necessary here, perform late for performance        
    tsincr = tsvel_matrix * span
    tscum = cumsum(tsincr, 2)
    # SB: convert zeros to nans not performed in matlab
    # SB: perform this after tsvel_matrix has been nan converted,
    # saves the step of comparing a large matrix (tsincr) to zero.
    # tscum = where(tscum == 0, nan, tscum) 

    return tsincr, tscum, tsvel_matrix


def time_series_by_rows(row, B0, SMFACTOR, SMORDER, ifg_data, mst, ncols,
                        nvelpar, PTHRESH, vcmt, TSMETHOD, INTERP):
    tsvel = np.empty(shape=(ncols, nvelpar), dtype=np.float32)
    for col in range(ncols):
        tsvel[col, :] = time_series_by_pixel(
            row, col, B0, SMFACTOR, SMORDER, ifg_data, mst, nvelpar,
            PTHRESH, vcmt, TSMETHOD, INTERP)

    return tsvel


def remove_rank_def_rows(B, nvelpar, ifgv, sel):
    q_var, r_var, e_var = qr(B, mode='economic', pivoting=True)
    licols = e_var[matrix_rank(B):nvelpar]
    [rmrow, rmcol] = where(B[:, licols] != 0)
    B = delete(B, rmrow, axis=0)
    ifgv = delete(ifgv, rmrow)
    sel = delete(sel, rmrow)
    return B, ifgv, sel, rmrow


def time_series_by_pixel(row, col, B0, SMFACTOR, SMORDER, ifg_data, mst,
                         nvelpar, PTHRESH, vcmt, method, interp):
    # check pixel for non-redundant ifgs
    sel = np.nonzero(mst[:, row, col])[0]  # trues in mst are chosen
    if len(sel) >= PTHRESH:
        ifgv = ifg_data[sel, row, col]
        # make design matrix, B
        B = B0[sel, :]
        if interp == 0:
            # remove rank deficient rows
            rmrow = asarray([0])  # dummy

            while len(rmrow) > 0:
                # the if else is introduced to prevent the following:
                # https://github.com/GeoscienceAustralia/PyRate/issues/69
                if B.shape[0] > 1:
                    B, ifgv, sel, rmrow = remove_rank_def_rows(
                        B, nvelpar, ifgv, sel)
                else:
                    return np.empty(nvelpar) * np.nan

            # Some epochs have been deleted; get valid epoch indices
            velflag = sum(abs(B), 0)
            # remove corresponding columns in design matrix
            B = B[:, ~np.isclose(velflag, 0.0)]
        else:
            velflag = np.ones(nvelpar)
        if method == 1:
            # Use Laplacian smoothing method
            tsvel = solve_ts_lap(nvelpar, velflag, ifgv, B,
                                 SMORDER, SMFACTOR, sel, vcmt)
        elif method == 2:
            # Use SVD method
            tsvel = solve_ts_svd(nvelpar, velflag, ifgv, B)
        else:
            raise ValueError("Unrecognised time series method")
        return tsvel
    else:
        return np.empty(nvelpar) * np.nan


def solve_ts_svd(nvelpar, velflag, ifgv, B):
    """
    Solve the linear least squares system using the SVD method.
    Similar to the SBAS method implemented by Berardino et al. 2002
    """
    # pre-allocate the velocity matrix
    tsvel = np.empty(nvelpar, dtype=np.float32) * np.nan
    # solve least squares equation using Moore-Penrose pseudoinverse
    tsvel[velflag != 0] = dot(pinv(B), ifgv)
    return tsvel


def solve_ts_lap(nvelpar, velflag, ifgv, B, SMORDER, SMFACTOR, sel, vcmt):
    """
    Solve the linear least squares system using the Finite Difference
    method using a Laplacian Smoothing operator.
    Similar to the method implemented by Schmidt and Burgmann 2003
    """
    # Laplacian observations number
    nlap = nvelpar - SMORDER  
    #  Laplacian smoothing coefficient matrix
    BLap0 = np.zeros(shape=(nlap, nvelpar))

    for i in range(nlap):
        if SMORDER == 1:
            BLap0[i, i:i+2] = [-1, 1]
        else:
            BLap0[i, i:i+3] = [1, -2, 1]

    # Scale the coefficients by Laplacian smoothing factor
    BLap0 *= SMFACTOR

    # Laplacian smoothing design matrix
    nvelleft = np.count_nonzero(velflag)
    nlap = nvelleft - SMORDER

    # constrain for the first and the last incremental
    BLap1 = -np.divide(np.ones(shape=nvelleft), nvelleft - 1)
    BLap1[0] = 1.0
    BLapn = -np.divide(np.ones(shape=nvelleft), nvelleft - 1)
    BLapn[-1] = 1.0

    BLap = np.empty(shape=(nlap + 2, nvelleft))
    BLap[0, :] = BLap1
    BLap[1:nlap + 1, :] = BLap0[0:nlap, 0:nvelleft]
    BLap[-1, :] = BLapn

    nlap += 2

    # add laplacian design matrix to existing design matrix
    B = np.concatenate((B, BLap), axis=0)

    # combine ifg and Laplacian smooth vector
    vLap = np.zeros(nlap)
    obsv = np.concatenate((ifgv, vLap), axis=0)

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
    wb = dot(w, B)
    wl = dot(w, obsv)
    x = dot(pinv(wb, rcond=1e-8), wl)

    # TODO: implement residuals and roughness calculations
    tsvel = np.empty(nvelpar, dtype=np.float32) * np.nan
    tsvel[~np.isclose(velflag, 0.0, atol=1e-8)] = x[:nvelleft]

    # TODO: implement uncertainty estimates (tserror) like in matlab code
    return tsvel


def plot_timeseries(tsincr, tscum, tsvel, output_dir):
    """
    Plot the results of the timeseries analysis.

    A very basic plotting function used in code development

    .. todo:: this should be moved and tidied up
    """

    nvelpar = len(tsincr[0, 0, :])
    for i in range(nvelpar):

        fig = plt.figure()
        imgplot = plt.imshow(tsincr[:, :, i])
        imgplot.set_clim(-10, 10)
        plt.colorbar(ticks=[-10, 0, 10], orientation ='horizontal')
        plt.draw()
        plt.savefig(output_dir + \
                    'tsincr_' + str(i) + '.png')
        plt.close()
        fig = None

        fig = plt.figure()
        imgplot = plt.imshow(tscum[:, :, i])
        imgplot.set_clim(-10, 10)
        plt.colorbar(ticks=[-10, 0, 10], orientation ='horizontal')
        plt.draw()
        plt.savefig(output_dir + \
                    'tscum_' + str(i) + '.png')
        plt.close()
        fig = None

        fig = plt.figure()
        imgplot = plt.imshow(tsvel[:, :, i])
        imgplot.set_clim(-50, 100)
        plt.colorbar(ticks=[-50, 50, 100], orientation ='horizontal')
        plt.draw()
        plt.savefig(output_dir + \
                    'tsvel_' + str(i) + '.png')
        plt.close()
        fig = None


def check_time_series_params(head, PTHRESH):
    """
    Validates time series parameter. head is any Ifg.
    """

    def missing_option_error(option):
        '''
        Internal convenience function for raising similar errors.
        '''
        msg = "Missing '%s' in configuration options" % option
        raise config.ConfigException(msg)

    if PTHRESH is None:
        missing_option_error(config.TIME_SERIES_PTHRESH)

    if PTHRESH < 0.0 or PTHRESH > 1000:
        raise ValueError(
            "minimum number of coherent observations for a pixel"
            " TIME_SERIES_PTHRESH setting must be >= 0.0 and <= 1000")



class TimeSeriesError(Exception):
    """
    Generic exception for time series errors.
    """

    pass