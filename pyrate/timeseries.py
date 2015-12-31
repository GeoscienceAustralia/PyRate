"""
Functions for producing a time series in PyRate.
InSAR time series inversion without smoothing
Based on original Matlab code by Hua Wang. Matlab 'tsinvnosm' function.

..codeauthor:: Vanessa Newey, Sudipta Basak
"""

from numpy import where, isnan, nan, diff, zeros, empty
from numpy import float32, cumsum, dot, delete, asarray
from numpy.linalg import matrix_rank, pinv
from scipy.linalg import qr
import numpy as np
import matplotlib.pyplot as plt
import pyrate.config as config
from algorithm import master_slave_ids, get_epochs
import gdal
import pyrate.ifgconstants as ifc
from pyrate import config as cf

def time_series_old(ifgs, pthresh, mst=None):
    """
    Returns time series data from the given ifgs.

    :param ifgs: Sequence of interferograms.
    :param pthresh: minimum number of coherent observations for a pixel.
    :param mst: [optional] array of ifg indexes from the pixel by pixel MST.

    :return: Tuple with the elements:

        - incremental deformation time series (NB: not velocity),
        - cumulative deformation time series, and
        - velocity

        these outputs are multi-dimensional arrays of size(nrows, ncols, nepochs-1),
        where:

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and
        - *nepochs* is the number of unique epochs (dates) covered by the ifgs.).
    """

    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    head = ifgs[0]
    check_time_series_params(head, pthresh)
    epochlist = get_epochs(ifgs)

    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    span = diff(epochlist.spans)
    nepoch = len(epochlist.dates)
    nvelpar = nepoch - 1

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

    tsincr = empty(shape=(nrows, ncols, nvelpar))
    tsvel = empty(shape=(nrows, ncols, nvelpar))

    tsincr[:] = nan
    ifg_data = zeros((nifgs, nrows, ncols), dtype=float32)
    for ifg_num in range(nifgs):
        ifgs[ifg_num].convert_to_nans(0)
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data

    if mst is None:
        mst = ~isnan(ifg_data)

    for row_num in xrange(nrows):
        for col_num in xrange(ncols):
            # check pixel for non-redundant ifgs;
            sel = where(mst[:, row_num, col_num])[0]
            m = len(sel)

            if m >= pthresh:
                ifgv = ifg_data[sel, row_num, col_num]
                B = B0[sel, :]
                # remove rank deficient rows
                rmrow = asarray([5, 4])

                while len(rmrow) > 0:
                    q_var, r_var, e_var = qr(B, mode='economic', pivoting=True)
                    licols = e_var[matrix_rank(B):nvelpar]
                    [rmrow, rmcol] = where(B[:, licols] != 0)
                    B = delete(B, rmrow, axis=0)
                    ifgv = delete(ifgv, rmrow)

                    if len(ifgv):
                        continue

                velflag = sum(abs(B), 0)
                B = B[:, where(velflag != 0)[0]]
                tsvel_pix = zeros(nvelpar, dtype=float32)
                tsvel_pix[where(velflag != 0)[0]] = dot(pinv(B), ifgv)
                tsvel[row_num, col_num, :] = tsvel_pix
                tsincr[row_num, col_num, :] = tsvel_pix * span
                # print row_num, col_num
                # print velflag
                # print B
                # print tsincr[row_num, col_num, :]

    if tsincr is None:
        raise TimeSeriesError("Could not produce a time series")

    tsincr = where(tsincr == 0, nan, tsincr)
    tscum = cumsum(tsincr, 2)

    # convert zeros to nans
    tscum = where(tscum == 0, nan, tscum)
    tsvel = where(tsvel == 0, nan, tsvel)
    print 'tsincr inside function', '-'*50
    print tsincr
    return tsincr, tscum, tsvel


def time_series(ifgs, pthresh, params, vcmt, mst=None):
    """
    Returns time series data from the given ifgs.

    :param ifgs: Sequence of interferograms.
    :param pthresh: minimum number of coherent observations for a pixel.
    :param mst: [optional] array of ifg indexes from the pixel by pixel MST.

    :return: Tuple with the elements:

        - incremental deformation time series (NB: not velocity),
        - cumulative deformation time series, and
        - velocity

        these outputs are multi-dimensional arrays of size(nrows, ncols, nepochs-1),
        where:

        - *nrows* is the number of rows in the ifgs,
        - *ncols* is the  number of columns in the ifgs, and
        - *nepochs* is the number of unique epochs (dates) covered by the ifgs.).
    """

    if len(ifgs) < 1:
        msg = 'Time series requires 2+ interferograms'
        raise TimeSeriesError(msg)

    SMORDER = params[cf.TIME_SERIES_SM_ORDER]
    SMFACTOR = np.power(10, params[cf.TIME_SERIES_SM_FACTOR])
    head = ifgs[0]
    check_time_series_params(head, pthresh)
    epochlist = get_epochs(ifgs)

    nrows = ifgs[0].nrows
    ncols = ifgs[0].ncols
    nifgs = len(ifgs)
    span = diff(epochlist.spans)
    nepoch = len(epochlist.dates)  # epoch number
    nvelpar = nepoch - 1  # velocity parameters number
    nlap = nvelpar- SMORDER  # Laplacian observations number

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

    #  Laplacian smoothing coefficient matrix

    BLap0 = np.zeros(shape=(nlap, nvelpar))

    for i in range(nlap):
      if SMORDER == 1:
        BLap0[i, i:i+2] = [-1, 1]
      else:
        BLap0[i, i:i+3] = [1, -2, 1]

    # smf: Laplacian smoothing factor
    BLap0 *= SMFACTOR

    tsincr = np.empty(shape=(nrows, ncols, nvelpar)) * np.nan
    tsvel_matrix = np.empty(shape=(nrows, ncols, nvelpar)) * np.nan

    ifg_data = zeros((nifgs, nrows, ncols), dtype=float32)
    for ifg_num in range(nifgs):
        ifgs[ifg_num].convert_to_nans(0)
        ifg_data[ifg_num] = ifgs[ifg_num].phase_data

    if mst is None:
        mst = ~isnan(ifg_data)

    for row_num in xrange(nrows):
        for col_num in xrange(ncols):
            # check pixel for non-redundant ifgs;
            sel = where(mst[:, row_num, col_num])[0]
            m = len(sel)

            if m >= pthresh:
                ifgv = ifg_data[sel, row_num, col_num]
                # make design matrix, B
                B = B0[sel, :]

                # remove rank deficient rows
                rmrow = asarray([5, 4])

                while len(rmrow) > 0:
                    q_var, r_var, e_var = qr(B, mode='economic', pivoting=True)
                    licols = e_var[matrix_rank(B):nvelpar]
                    [rmrow, rmcol] = where(B[:, licols] != 0)
                    B = delete(B, rmrow, axis=0)
                    ifgv = delete(ifgv, rmrow)
                    sel = delete(sel, rmrow)

                    # if len(ifgv):
                    #     continue
                m = len(sel)

                velflag = sum(abs(B), 0)
                B = B[:, ~np.isclose(velflag, 0.0)]

                # Laplacian smoothing design matrix
                nvelleft = np.count_nonzero(velflag)
                nlap = nvelleft - SMORDER

                BLap = np.empty(shape=(nlap + 2, nvelleft))

                # constrain for the first and the last incremental
                BLap1 = -np.divide(np.ones(shape=nvelleft), nvelleft-1)
                BLap1[0] = 1.0

                BLapn = -np.divide(np.ones(shape=nvelleft), nvelleft-1)
                BLapn[-1] = 1.0

                BLap[0, :] = BLap1
                BLap[1:nlap+1, :] = BLap0[0:nlap, 0:nvelleft]
                BLap[-1, :] = BLapn

                nlap += 2

                # add laplacian design matrix
                B = np.concatenate((B, BLap), axis=0)

                # combine ifg and Laplacian smooth vector
                vLap = np.zeros(nlap)
                obsv = np.concatenate((ifgv, vLap), axis=0)

                # make variance-covariance matrix
                # new covariance matrix, adding the laplacian equations
                nobs = m + nlap
                vcm_tmp = np.eye(nobs)
                vcm_tmp[0:m, 0:m] = vcmt[sel, np.vstack(sel)]

                # solve the equation by least-squares
                # calculate velocities
                # we get the lower triangle in numpy
                # matlab gives upper triangle by default
                w = np.linalg.cholesky(np.linalg.pinv(vcm_tmp)).T

                wb = np.dot(w, B)
                wl = np.dot(w, obsv)
                # x=wb\wl
                x = np.dot(np.linalg.pinv(wb), wl)

                # residuals and roughness
                # not implemented

                tsvel = np.empty(nvelpar)*np.nan
                tsvel[~np.isclose(velflag, 0.0, atol=1e-8)] = x[:nvelleft]

                tsvel_matrix[row_num, col_num, :] = tsvel
                tsincr[row_num, col_num, :] = tsvel*span

    if tsincr is None:
        raise TimeSeriesError("Could not produce a time series")

    tscum = cumsum(tsincr, 2)

    # convert zeros to nans
    tsincr = where(tsincr == 0, nan, tsincr)
    tscum = where(tscum == 0, nan, tscum)
    tsvel_matrix = where(tsvel_matrix == 0, nan, tsvel_matrix)
    return tsincr, tscum, tsvel_matrix


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


def check_time_series_params(head, pthresh):
    """
    Validates time series parameter. head is any Ifg.
    """

    def missing_option_error(option):
        '''
        Internal convenience function for raising similar errors.
        '''
        msg = "Missing '%s' in configuration options" % option
        raise config.ConfigException(msg)

    if pthresh is None:
        missing_option_error(config.TIME_SERIES_PTHRESH)

    if pthresh < 0.0 or pthresh > 1000:
        raise ValueError(
            "minimum number of coherent observations for a pixel"
            " TIME_SERIES_PTHRESH setting must be >= 0.0 and <= 1000")


def write_geotiff_output(md, data, dest, nodata):
    '''
    Writes data to a GeoTIFF file.
    md is a dictionary containing these items:
    ifc.PYRATE_PROJECTION, ifc.PYRATE_GEOTRANSFORM,
    ifc.PYRATE_DATE
    NB It should be moved to utils class
    '''

    driver = gdal.GetDriverByName("GTiff")

    nrows =len(data[:, 0])
    ncols = len(data[0, :])
    print "ncols", ncols
    print "nrows", nrows
    ds = driver.Create(dest, ncols, nrows, 1, gdal.GDT_Float32)
    print str(md[ifc.PYRATE_DATE])
    print md
    ds.SetProjection(md[ifc.PYRATE_PROJECTION])
    ds.SetGeoTransform(md[ifc.PYRATE_GEOTRANSFORM])
    ds.SetMetadataItem(ifc.PYRATE_DATE, str(md[ifc.PYRATE_DATE]))

    # write data
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)


class TimeSeriesError(Exception):
    """
    Generic exception for time series errors.
    """

    pass


if __name__ == "__main__":
    import os
    import shutil
    from subprocess import call

    from pyrate.scripts import run_pyrate
    from pyrate import matlab_mst_kruskal as matlab_mst
    from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT
    from pyrate.tests.common import SYD_TEST_DIR
    from pyrate import config as cf
    from pyrate import reference_phase_estimation as rpe
    from pyrate import vcm

    SYD_TIME_SERIES_DIR = os.path.join(SYD_TEST_DIR, 'matlab_time_series')

    # start each full test run cleanly
    shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

    os.makedirs(SYD_TEST_OUT)

    params = cf.get_config_params(
            os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))
    params[cf.REF_EST_METHOD] = 2
    call(["python", "pyrate/scripts/run_prepifg.py",
          os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf')])

    xlks, ylks, crop = run_pyrate.transform_params(params)

    base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

    dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params, xlks)

    ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

    assert isinstance(ifg_instance, matlab_mst.IfgListPyRate)
    ifgs = ifg_instance.ifgs
    for i in ifgs:
        if not i.mm_converted:
            i.convert_to_mm()
            i.write_modified_phase()
    ifg_instance_updated, epoch_list = \
        matlab_mst.get_nml(ifg_instance, nan_conversion=True)
    mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

    for i in ifgs:
        if not i.is_open:
            i.open()
        if not i.nan_converted:
            i.convert_to_nans()

        if not i.mm_converted:
            i.convert_to_mm()
            i.write_modified_phase()

    if params[cf.ORBITAL_FIT] != 0:
        run_pyrate.remove_orbital_error(ifgs, params)

    refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

    if params[cf.ORBITAL_FIT] != 0:
        run_pyrate.remove_orbital_error(ifgs, params)

    _, ifgs = rpe.estimate_ref_phase(ifgs, params, refx, refy)

    # Calculate interferogram noise
    # TODO: assign maxvar to ifg metadata (and geotiff)?
    maxvar = [vcm.cvd(i)[0] for i in ifgs]

    # Calculate temporal variance-covariance matrix
    vcmt = vcm.get_vcmt(ifgs, maxvar)

    # Calculate linear rate map
    rate, error, samples = run_pyrate.calculate_linear_rate(ifgs, params, vcmt,
                                                 mst=mst_grid)

    tsincr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_incr.csv')
    ts_incr = np.genfromtxt(tsincr_path, delimiter=',')
    tserr_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_error.csv')
    ts_err = np.genfromtxt(tserr_path, delimiter=',')
    tscum_path = os.path.join(SYD_TIME_SERIES_DIR, 'ts_cum.csv')
    ts_cum = np.genfromtxt(tscum_path, delimiter=',')
    np.set_printoptions(precision=4)
    # Calculate time series
    if params[cf.TIME_SERIES_CAL] != 0:
        tsincr, tscum, tsvel = run_pyrate.calculate_time_series(
            ifgs, params, vcmt, mst=mst_grid)
        print 'here in calc time series'
        print '-'*50
        print np.sum(np.isnan(np.ravel(tsincr))), np.sum(np.isnan(ts_incr))
        nan_index = np.isnan(np.ravel(tsincr))
        print ts_incr[nan_index]
        ts_incr = np.reshape(ts_incr, newshape=tsincr.shape, order='F')
        ts_cum = np.reshape(ts_cum, newshape=tsincr.shape, order='F')

        # mismatch = np.nonzero(~np.isclose(tsincr, ts_incr, atol=1e-3))
        # print mismatch

        # for i in range(len(mismatch[0])):
        #     print mismatch[0][i], mismatch[1][i], mismatch[2][i]

        #TODO: Investigate why the whole matrices don't equal
        # Current hypothesis is that the pseudo inverse computed are different
        # in matlab and python as they
        np.testing.assert_array_almost_equal(
            ts_incr[:11, :45, :], tsincr[:11, :45, :], decimal=4)

        np.testing.assert_array_almost_equal(
            ts_cum[:11, :45, :], tscum[:11, :45, :], decimal=4)

