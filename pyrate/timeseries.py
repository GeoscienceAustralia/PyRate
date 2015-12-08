"""
Functions for producing a time series in PyRate.
InSAR time series inversion without smoothing
Based on original Matlab code by Hua Wang. Matlab 'tsinvnosm' function.

..codeauthor:: Vanessa Newey
"""

from numpy import where, isnan, nan, diff, zeros, empty
from numpy import float32, cumsum, dot, delete, asarray
from numpy.linalg import matrix_rank, pinv
from scipy.linalg import qr
import matplotlib.pyplot as plt
import pyrate.config as config
from algorithm import master_slave_ids, get_epochs
import gdal
import pyrate.ifgconstants as ifc


def time_series(ifgs, pthresh, mst=None):
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
#               remove rank deficient rows
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
                B = B[:, where(velflag!=0)[0]]
                tsvel_pix = zeros(nvelpar, dtype=float32)
                tsvel_pix[where(velflag!=0)[0]] = dot(pinv(B), ifgv)
                tsvel[row_num, col_num, :] = tsvel_pix
                tsincr[row_num, col_num, :] = tsvel_pix * span

    if tsincr is None:
        raise TimeSeriesError("Could not produce a time series")

    tsincr = where(tsincr == 0, nan, tsincr)
    tscum = cumsum(tsincr, 2)

    #convert zeros to nans
    tscum = where(tscum == 0, nan, tscum)
    tsvel = where(tsvel == 0, nan, tsvel)
    return tsincr, tscum, tsvel


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
        raise ValueError("minimum number of coherent observations for a pixel" +
                          + " TIME_SERIES_PTHRESH setting must be >= 0.0 and <= 1000")


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
