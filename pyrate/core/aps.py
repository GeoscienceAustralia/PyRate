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
This Python module implements a spatio-temporal filter method
for correcting interferograms for atmospheric phase screen (APS)
signals.
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import os
from copy import deepcopy
from collections import OrderedDict
from typing import List
import numpy as np
from numpy import isnan
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift
from scipy.interpolate import griddata
from pyrate.core.logger import pyratelogger as log

from pyrate.core import shared, ifgconstants as ifc, mpiops, config as cf
from pyrate.core.covariance import cvd_from_phase, RDist
from pyrate.core.algorithm import get_epochs
from pyrate.core.shared import Ifg
from pyrate.core.timeseries import time_series
from pyrate.merge import assemble_tiles
from pyrate.configuration import MultiplePaths, Configuration


def wrap_spatio_temporal_filter(params):
    """
    A wrapper for the spatio-temporal filter so it can be tested.
    See docstring for spatio_temporal_filter.
    """
    if params[cf.APSEST]:
        log.info('Doing APS spatio-temporal filtering')
    else:
        log.info('APS spatio-temporal filtering not required')
        return
    tiles = params[cf.TILES]
    preread_ifgs = params[cf.PREREAD_IFGS]
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]

    # perform some checks on existing ifgs
    log.debug('Checking APS correction status')
    if mpiops.run_once(shared.check_correction_status, ifg_paths, ifc.PYRATE_APS_ERROR):
        log.debug('Finished APS correction')
        return  # return if True condition returned

    aps_error_files_on_disc = [MultiplePaths.aps_error_path(i, params) for i in ifg_paths]
    if all(a.exists() for a in aps_error_files_on_disc):
        log.warning("Reusing APS errors from previous run!!!")
        for ifg_path, a in mpiops.array_split(list(zip(ifg_paths, aps_error_files_on_disc))):
            phase = np.load(a)
            _save_aps_corrected_phase(ifg_path, phase)
    else:
        tsincr = _calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles)
        mpiops.comm.barrier()

        spatio_temporal_filter(tsincr, ifg_paths, params, preread_ifgs)
    mpiops.comm.barrier()
    shared.save_numpy_phase(ifg_paths, params)


def spatio_temporal_filter(tsincr, ifg_paths, params, preread_ifgs):
    """
    Applies a spatio-temporal filter to remove the atmospheric phase screen
    (APS) and saves the corrected interferograms. Before performing this step,
    the time series is computed using the SVD method. This function then
    performs temporal and spatial filtering.

    :param ndarray tsincr: incremental time series array of size
                (ifg.shape, nepochs-1)
    :param list ifg_paths: List of interferogram file paths
    :param dict params: Dictionary of configuration parameter
    :param dict preread_ifgs: Dictionary of shared.PrereadIfg class instances

    :return: None, corrected interferograms are saved to disk
    """
    ifg = Ifg(ifg_paths[0])  # just grab any for parameters in slpfilter
    ifg.open()
    epochlist = mpiops.run_once(get_epochs, preread_ifgs)[0]
    ts_hp = temporal_high_pass_filter(tsincr, epochlist, params)
    ts_aps = spatial_low_pass_filter(ts_hp, ifg, params)
    tsincr -= ts_aps

    _ts_to_ifgs(tsincr, preread_ifgs, params)
    ifg.close()


def _calc_svd_time_series(ifg_paths, params, preread_ifgs,
                          tiles: List[shared.Tile]):
    """
    Helper function to obtain time series for spatio-temporal filter
    using SVD method
    """
    # Is there other existing functions that can perform this same job?
    log.info('Calculating incremental time series via SVD method for APS '
             'correction')
    # copy params temporarily
    new_params = deepcopy(params)
    new_params[cf.TIME_SERIES_METHOD] = 2  # use SVD method

    process_tiles = mpiops.array_split(tiles)

    nvels = None
    for t in process_tiles:
        log.debug(f'Calculating time series for tile {t.index} during APS '
                  f'correction')
        ifgp = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
        mst_tile = np.load(Configuration.mst_path(params, t.index))
        tsincr = time_series(ifgp, new_params, vcmt=None, mst=mst_tile)[0]
        np.save(file=os.path.join(params[cf.TMPDIR],
                f'tsincr_aps_{t.index}.npy'), arr=tsincr)
        nvels = tsincr.shape[2]

    nvels = mpiops.comm.bcast(nvels, root=0)
    mpiops.comm.barrier()
    # need to assemble tsincr from all processes
    tsincr_g = _assemble_tsincr(ifg_paths, params, preread_ifgs, tiles, nvels)
    log.debug('Finished calculating time series for spatio-temporal filter')
    return tsincr_g


def _assemble_tsincr(ifg_paths, params, preread_ifgs, tiles, nvels):
    """
    Helper function to reconstruct time series images from tiles
    """
    # pre-allocate dest 3D array
    shape = preread_ifgs[ifg_paths[0]].shape
    tsincr_p = {}
    process_nvels = mpiops.array_split(range(nvels))
    for i in process_nvels:
        tsincr_p[i] = assemble_tiles(shape, params[cf.TMPDIR], tiles,
                                     out_type='tsincr_aps', index=i)
    tsincr_g = shared.join_dicts(mpiops.comm.allgather(tsincr_p))
    return np.dstack([v[1] for v in sorted(tsincr_g.items())])


def _ts_to_ifgs(tsincr, preread_ifgs, params):
    """
    Function that converts an incremental displacement time series into
    interferometric phase observations. Used to re-construct an interferogram
    network from a time series.

    :param ndarray tsincr: incremental time series array of size
                (ifg.shape, nepochs-1)
    :param dict preread_ifgs: Dictionary of shared.PrereadIfg class instances

    :return: None, interferograms are saved to disk
    """
    log.debug('Reconstructing interferometric observations from time series')
    ifgs = list(OrderedDict(sorted(preread_ifgs.items())).values())
    _, n = mpiops.run_once(get_epochs, ifgs)
    index_first, index_second = n[:len(ifgs)], n[len(ifgs):]

    num_ifgs_tuples = mpiops.array_split(list(enumerate(ifgs)))
    num_ifgs_tuples = [(int(num), ifg) for num, ifg in num_ifgs_tuples]

    for i, ifg in num_ifgs_tuples:
        aps_error_on_disc = MultiplePaths.aps_error_path(ifg.tmp_path, params)
        phase = np.sum(tsincr[:, :, index_first[i]: index_second[i]], axis=2)
        np.save(file=aps_error_on_disc, arr=phase)
        _save_aps_corrected_phase(ifg.tmp_path, phase)


def _save_aps_corrected_phase(ifg_path, phase):
    """
    Save (update) interferogram metadata and phase data after
    spatio-temporal filter (APS) correction.
    """
    ifg = Ifg(ifg_path)
    ifg.open(readonly=False)
    ifg.phase_data[~np.isnan(ifg.phase_data)] = phase[~np.isnan(ifg.phase_data)]
    # set aps tags after aps error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_APS_ERROR, ifc.APS_REMOVED)
    ifg.write_modified_phase()
    ifg.close()


def spatial_low_pass_filter(ts_hp, ifg, params):
    """
    Filter time series data spatially using a Gaussian low-pass
    filter defined by a cut-off distance. If the cut-off distance is
    defined as zero in the parameters dictionary then it is calculated for
    each time step using the pyrate.covariance.cvd_from_phase method.

    :param ndarray ts_hp: Array of time series data, the result of a temporal
                       high-pass filter operation. shape (ifg.shape, n_epochs)
    :param shared.Ifg instance ifg: interferogram object
    :param dict params: Dictionary of configuration parameters

    :return: ts_lp: low-pass filtered time series data of shape
                    (ifg.shape, n_epochs)
    :rtype: ndarray
    """
    log.info('Applying spatial low-pass filter')
    if params[cf.SLPF_NANFILL] == 0:
        ts_hp[np.isnan(ts_hp)] = 0  # need it here for cvd and fft
    else:
        # optionally interpolate, operation is inplace
        _interpolate_nans(ts_hp, params[cf.SLPF_NANFILL_METHOD])

    nvels = ts_hp.shape[2]
    cutoff = params[cf.SLPF_CUTOFF]
    if cutoff == 0:
        r_dist = RDist(ifg)() # only needed for cvd_for_phase
    else:
        r_dist = None
        log.info(f'Gaussian spatial filter cutoff is {cutoff:.3f} km for all '
                 f'{nvels} time-series images')

    process_nvel = mpiops.array_split(range(nvels))
    process_ts_lp = {}

    for i in process_nvel:
        process_ts_lp[i] = _slpfilter(ts_hp[:, :, i], ifg, r_dist, cutoff)

    ts_lp_d = shared.join_dicts(mpiops.comm.allgather(process_ts_lp))
    ts_lp = np.dstack([v[1] for v in sorted(ts_lp_d.items())])
    log.debug('Finished applying spatial low pass filter')
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
        method=method)
    a[np.isnan(a)] = 0  # zero fill boundary/edge nans


def _slpfilter(phase, ifg, r_dist, cutoff):
    """
    Wrapper function for spatial low pass filter
    """
    if np.all(np.isnan(phase)):  # return for nan matrix
        return phase

    if cutoff == 0:
        _, alpha = cvd_from_phase(phase, ifg, r_dist, calc_alpha=True)
        cutoff = 1.0/alpha
        log.info(f'Gaussian spatial filter cutoff is {cutoff:.3f} km')

    return gaussian_spatial_filter(phase, cutoff, ifg.x_size, ifg.y_size)


def gaussian_spatial_filter(image, cutoff, x_size, y_size):
    """
    Function to apply a Gaussian spatial low-pass filter to a 2D image with
    unequal pixel resolution in x and y dimensions.

    :param ndarray image: 2D image to be filtered
    :param float cutoff: filter cutoff in kilometres
    :param float x_size: pixel size in x dimension, in metres
    :param float y_size: pixel size in y dimension, in metres

    :return: out: Gaussian low-pass filtered 2D image
    :rtype: ndarray
    """
    rows, cols = image.shape
    cx = np.floor(cols/2)
    cy = np.floor(rows/2)
    # fft for the input image
    imf = fftshift(fft2(image))
    # calculate distance
    [xx, yy] = np.meshgrid(range(cols), range(rows))
    xx = (xx - cx) * x_size  # these are in meters as x_size in metres
    yy = (yy - cy) * y_size
    dist = np.sqrt(xx ** 2 + yy ** 2)/ ifc.METRE_PER_KM # change m to km

    # Apply Gaussian smoothing kernel
    H = _kernel(dist, cutoff)
    outf = imf * H
    out = np.real(ifft2(ifftshift(outf)))
    out[np.isnan(image)] = np.nan # re-apply nans to output image
    return out


# TODO: use tiles here and distribute amongst processes
def temporal_high_pass_filter(tsincr, epochlist, params):
    """
    Isolate high-frequency components of time series data by subtracting
    low-pass components obtained using a Gaussian filter defined by a
    cut-off time period (in days).

    :param ndarray tsincr: Array of incremental time series data of shape
                (ifg.shape, n_epochs)
    :param list epochlist: List of shared.EpochList class instances
    :param dict params: Dictionary of configuration parameters

    :return: ts_hp: filtered high frequency time series data,
                    shape (ifg.shape, nepochs)
    :rtype: ndarray
    """
    log.info('Applying temporal high-pass filter')
    threshold = params[cf.TLPF_PTHR]
    cutoff_day = params[cf.TLPF_CUTOFF]
    if cutoff_day < 1 or type(cutoff_day) != int:
        raise ValueError(f'tlpf_cutoff must be an integer greater than or '
                         f'equal to 1 day. Value provided = {cutoff_day}')

    # convert cutoff in days to years
    cutoff_yr = cutoff_day / ifc.DAYS_PER_YEAR
    log.info(f'Gaussian temporal filter cutoff is {cutoff_day} days '
             f'({cutoff_yr:.4f} years)')

    intv = np.diff(epochlist.spans)  # time interval for the neighboring epochs
    span = epochlist.spans[: tsincr.shape[2]] + intv/2  # accumulated time
    rows, cols = tsincr.shape[:2]

    tsfilt_row = {}
    process_rows = mpiops.array_split(list(range(rows)))

    for r in process_rows:
        tsfilt_row[r] = np.empty(tsincr.shape[1:], dtype=np.float32) * np.nan
        for j in range(cols):
            # Result of gaussian filter is low frequency time series
            tsfilt_row[r][j, :] = gaussian_temporal_filter(tsincr[r, j, :],
                                                     cutoff_yr, span, threshold)

    tsfilt_combined = shared.join_dicts(mpiops.comm.allgather(tsfilt_row))
    tsfilt = np.array([v[1] for v in tsfilt_combined.items()])
    log.debug("Finished applying temporal high-pass filter")
    # Return the high-pass time series by subtracting low-pass result from input
    return tsincr - tsfilt


def gaussian_temporal_filter(tsincr, cutoff, span, thr):
    """
    Function to apply a Gaussian temporal low-pass filter to a 1D time-series
    vector for one pixel with irregular temporal sampling.

    :param ndarray tsincr: 1D time-series vector to be filtered
    :param float cutoff: filter cutoff in years
    :param ndarray span: 1D vector of cumulative time spans, in years
    :param int thr: threshold for non-nan values in tsincr

    :return: ts_lp: low-pass filtered time series vector
    :rtype: ndarray
    """
    nanmat = ~isnan(tsincr)
    sel = np.nonzero(nanmat)[0]  # don't select if nan
    ts_lp = np.empty(tsincr.shape, dtype=np.float32) * np.nan
    m = len(sel)
    if m >= thr:
        for k in range(m):
            yr = span[sel] - span[sel[k]]
            # apply Gaussian smoothing kernel
            wgt = _kernel(yr, cutoff)
            wgt /= np.sum(wgt)
            ts_lp[sel[k]] = np.sum(tsincr[sel] * wgt)

    return ts_lp

def _kernel(x, cutoff):
    """
    Gaussian low-pass filter kernel
    """
    return np.exp(-0.5 * (x / cutoff) ** 2)
