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

import pyrate.constants as C
from pyrate.core.logger import pyratelogger as log

from pyrate.core import shared, ifgconstants as ifc, mpiops
from pyrate.core.covariance import cvd_from_phase, RDist
from pyrate.core.algorithm import get_epochs
from pyrate.core.shared import Ifg, Tile, EpochList, nan_and_mm_convert
from pyrate.core.timeseries import time_series
from pyrate.merge import assemble_tiles
from pyrate.configuration import MultiplePaths, Configuration


def spatio_temporal_filter(params: dict) -> None:
    """
    Applies a spatio-temporal filter to remove the atmospheric phase screen
    (APS) and saves the corrected interferograms. Firstly the incremental
    time series is computed using the SVD method, before a cascade of temporal
    then spatial Gaussian filters is applied. The resulting APS corrections are
    saved to disc before being subtracted from each interferogram.

    :param params: Dictionary of PyRate configuration parameters.
    """
    if params[C.APSEST]:
        log.info('Doing APS spatio-temporal filtering')
    else:
        log.info('APS spatio-temporal filtering not required')
        return
    tiles = params[C.TILES]
    preread_ifgs = params[C.PREREAD_IFGS]
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]

    # perform some checks on existing ifgs
    log.debug('Checking APS correction status')
    if mpiops.run_once(shared.check_correction_status, ifg_paths, ifc.PYRATE_APS_ERROR):
        log.debug('Finished APS correction')
        return  # return if True condition returned

    aps_paths = [MultiplePaths.aps_error_path(i, params) for i in ifg_paths]
    if all(a.exists() for a in aps_paths):
        log.warning('Reusing APS errors from previous run')
        _apply_aps_correction(ifg_paths, aps_paths, params)
        return

    # obtain the incremental time series using SVD
    tsincr = _calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles)
    mpiops.comm.barrier()

    # get lists of epochs and ifgs
    ifgs = list(OrderedDict(sorted(preread_ifgs.items())).values())
    epochlist = mpiops.run_once(get_epochs, ifgs)[0]

    # first perform temporal high pass filter
    ts_hp = temporal_high_pass_filter(tsincr, epochlist, params)

    # second perform spatial low pass filter to obtain APS correction in ts domain
    ifg = Ifg(ifg_paths[0])  # just grab any for parameters in slpfilter
    ifg.open()
    ts_aps = spatial_low_pass_filter(ts_hp, ifg, params)
    ifg.close()

    # construct APS corrections for each ifg
    _make_aps_corrections(ts_aps, ifgs, params)

    # apply correction to ifgs and save ifgs to disc.
    _apply_aps_correction(ifg_paths, aps_paths, params)

    # update/save the phase_data in the tiled numpy files
    shared.save_numpy_phase(ifg_paths, params)


def _calc_svd_time_series(ifg_paths: List[str], params: dict, preread_ifgs: dict,
                          tiles: List[Tile]) -> np.ndarray:
    """
    Helper function to obtain time series for spatio-temporal filter
    using SVD method
    """
    # Is there other existing functions that can perform this same job?
    log.info('Calculating incremental time series via SVD method for APS '
             'correction')
    # copy params temporarily
    new_params = deepcopy(params)
    new_params[C.TIME_SERIES_METHOD] = 2  # use SVD method

    process_tiles = mpiops.array_split(tiles)

    nvels = None
    for t in process_tiles:
        log.debug(f'Calculating time series for tile {t.index} during APS '
                  f'correction')
        ifgp = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
        mst_tile = np.load(Configuration.mst_path(params, t.index))
        tsincr = time_series(ifgp, new_params, vcmt=None, mst=mst_tile)[0]
        np.save(file=os.path.join(params[C.TMPDIR],
                                  f'tsincr_aps_{t.index}.npy'), arr=tsincr)
        nvels = tsincr.shape[2]

    nvels = mpiops.comm.bcast(nvels, root=0)
    mpiops.comm.barrier()
    # need to assemble tsincr from all processes
    tsincr_g = _assemble_tsincr(ifg_paths, params, preread_ifgs, tiles, nvels)
    log.debug('Finished calculating time series for spatio-temporal filter')
    return tsincr_g


def _assemble_tsincr(ifg_paths: List[str], params: dict, preread_ifgs: dict,
                     tiles: List[Tile], nvels: np.float32) -> np.ndarray:
    """
    Helper function to reconstruct time series images from tiles
    """
    # pre-allocate dest 3D array
    shape = preread_ifgs[ifg_paths[0]].shape
    tsincr_p = {}
    process_nvels = mpiops.array_split(range(nvels))
    for i in process_nvels:
        tsincr_p[i] = assemble_tiles(shape, params[C.TMPDIR], tiles,
                                     out_type='tsincr_aps', index=i)
    tsincr_g = shared.join_dicts(mpiops.comm.allgather(tsincr_p))
    return np.dstack([v[1] for v in sorted(tsincr_g.items())])


def _make_aps_corrections(ts_aps: np.ndarray, ifgs: List[Ifg], params: dict) -> None:
    """
    Function to convert the time series APS filter output into interferometric
    phase corrections and save them to disc.

    :param ts_aps: Incremental APS time series array.
    :param ifgs:   List of Ifg class objects.
    :param params: Dictionary of PyRate configuration parameters.
    """
    log.debug('Reconstructing interferometric observations from time series')
    # get first and second image indices 
    _ , n = mpiops.run_once(get_epochs, ifgs)
    index_first, index_second = n[:len(ifgs)], n[len(ifgs):]

    num_ifgs_tuples = mpiops.array_split(list(enumerate(ifgs)))
    for i, ifg in [(int(num), ifg) for num, ifg in num_ifgs_tuples]:
        # sum time slice data from first to second epoch
        ifg_aps = np.sum(ts_aps[:, :, index_first[i]: index_second[i]], axis=2)
        aps_error_on_disc = MultiplePaths.aps_error_path(ifg.tmp_path, params)
        np.save(file=aps_error_on_disc, arr=ifg_aps) # save APS as numpy array

    mpiops.comm.barrier()


def _apply_aps_correction(ifg_paths: List[str], aps_paths: List[str], params: dict) -> None:
    """
    Function to read and apply (subtract) APS corrections from interferogram data.
    """
    for ifg_path, aps_path in mpiops.array_split(list(zip(ifg_paths, aps_paths))):
        # read the APS correction from numpy array
        aps_corr = np.load(aps_path)
        # open the Ifg object
        ifg = Ifg(ifg_path)
        ifg.open(readonly=False)
        # convert NaNs and convert to mm
        nan_and_mm_convert(ifg, params)
        # subtract the correction from the ifg phase data
        ifg.phase_data[~np.isnan(ifg.phase_data)] -= aps_corr[~np.isnan(ifg.phase_data)]
        # set meta-data tags after aps error correction
        ifg.dataset.SetMetadataItem(ifc.PYRATE_APS_ERROR, ifc.APS_REMOVED)
        # write phase data to disc and close ifg.
        ifg.write_modified_phase()
        ifg.close()


def spatial_low_pass_filter(ts_hp: np.ndarray, ifg: Ifg, params: dict) -> np.ndarray:
    """
    Filter time series data spatially using a Gaussian low-pass
    filter defined by a cut-off distance. If the cut-off distance is
    defined as zero in the parameters dictionary then it is calculated for
    each time step using the pyrate.covariance.cvd_from_phase method.
    :param ts_hp: Array of temporal high-pass time series data, shape (ifg.shape, n_epochs)
    :param ifg: pyrate.core.shared.Ifg Class object.
    :param params: Dictionary of PyRate configuration parameters.
    :return: ts_lp: Low-pass filtered time series data of shape (ifg.shape, n_epochs).
    """
    log.info('Applying spatial low-pass filter')

    nvels = ts_hp.shape[2]
    cutoff = params[C.SLPF_CUTOFF]
    # nanfill = params[cf.SLPF_NANFILL]
    # fillmethod = params[cf.SLPF_NANFILL_METHOD]
    if cutoff == 0:
        r_dist = RDist(ifg)()  # only needed for cvd_for_phase
    else:
        r_dist = None
        log.info(f'Gaussian spatial filter cutoff is {cutoff:.3f} km for all '
                 f'{nvels} time-series images')

    process_nvel = mpiops.array_split(range(nvels))
    process_ts_lp = {}

    for i in process_nvel:
        process_ts_lp[i] = _slpfilter(ts_hp[:, :, i], ifg, r_dist, params)

    ts_lp_d = shared.join_dicts(mpiops.comm.allgather(process_ts_lp))
    ts_lp = np.dstack([v[1] for v in sorted(ts_lp_d.items())])
    log.debug('Finished applying spatial low pass filter')
    return ts_lp


def _interpolate_nans_2d(arr: np.ndarray, method: str) -> None:
    """
    In-place array interpolation and NaN-fill using scipy.interpolation.griddata.
    :param arr: 2D ndarray to be interpolated.
    :param method: Method; one of 'nearest', 'linear', and 'cubic'.
    """
    log.debug(f'Interpolating array with "{method}" method')
    r, c = np.indices(arr.shape)
    arr[np.isnan(arr)] = griddata(
        (r[~np.isnan(arr)], c[~np.isnan(arr)]),  # points we know
        arr[~np.isnan(arr)],  # values we know
        (r[np.isnan(arr)], c[np.isnan(arr)]),  # points to interpolate
        method=method, fill_value=0)


def _slpfilter(phase: np.ndarray, ifg: Ifg, r_dist: float, params: dict) -> np.ndarray:
    """
    Wrapper function for spatial low pass filter
    """
    cutoff = params[C.SLPF_CUTOFF]
    nanfill = params[C.SLPF_NANFILL]
    fillmethod = params[C.SLPF_NANFILL_METHOD]

    if np.all(np.isnan(phase)):  # return for nan matrix
        return phase

    if cutoff == 0:
        _, alpha = cvd_from_phase(phase, ifg, r_dist, calc_alpha=True)
        cutoff = 1.0 / alpha
        log.info(f'Gaussian spatial filter cutoff is {cutoff:.3f} km')

    return gaussian_spatial_filter(phase, cutoff, ifg.x_size, ifg.y_size, nanfill, fillmethod)


def gaussian_spatial_filter(image: np.ndarray, cutoff: float, x_size: float,
                            y_size: float, nanfill: bool = True,
                            fillmethod: str = 'nearest') -> np.ndarray:
    """
    Function to apply a Gaussian spatial low-pass filter to a 2D image with
    unequal pixel resolution in x and y dimensions. Performs filtering in the
    Fourier domain. Any NaNs in the image are interpolated prior to Fourier
    transformation, with NaNs being replaced in to the filtered output image.
    :param image: 2D image to be filtered
    :param cutoff: filter cutoff in kilometres
    :param x_size: pixel size in x dimension, in metres
    :param y_size: pixel size in y dimension, in metres
    :param nanfill: interpolate image to fill NaNs
    :param fillmethod: interpolation method ('nearest', 'cubic', or 'linear')
    :return: filt: Gaussian low-pass filtered 2D image
    """
    # create NaN mask of image
    mask = np.isnan(image)
    # in-place nearest-neighbour interpolation to fill NaNs
    # nearest neighbour will fill values outside the convex hull
    if nanfill:
        _interpolate_nans_2d(image, fillmethod)

    rows, cols = image.shape
    pad = 4096
    # pad the image to a large square array.
    # TODO: implement variable padding dependent on image size
    im = np.pad(image, ((0, pad - rows), (0, pad - cols)), 'constant')
    # fast fourier transform of the input image
    imf = fftshift(fft2(im))

    # calculate centre coords of image
    cx = np.floor(pad / 2)
    cy = np.floor(pad / 2)
    # calculate distance array
    [xx, yy] = np.meshgrid(range(pad), range(pad))
    xx = (xx - cx) * x_size  # these are in meters as x_size in metres
    yy = (yy - cy) * y_size
    dist = np.sqrt(xx ** 2 + yy ** 2) / ifc.METRE_PER_KM  # change m to km

    # Estimate sigma value for Gaussian kernel function in spectral domain
    # by converting cutoff distance to wavenumber and applying a scaling
    # factor based on fixed kernel window size. 
    sigma = np.std(dist) * (1 / cutoff)
    # Calculate kernel weights
    wgt = _kernel(dist, sigma)
    # Apply Gaussian smoothing kernel
    outf = imf * wgt
    # Inverse Fourier transform
    out = np.real(ifft2(ifftshift(outf)))
    filt = out[:rows, :cols]  # grab non-padded part
    filt[mask] = np.nan  # re-insert nans in output image
    return filt


# TODO: use tiles here and distribute amongst processes
def temporal_high_pass_filter(tsincr: np.ndarray, epochlist: EpochList,
                              params: dict) -> np.ndarray:
    """
    Isolate high-frequency components of time series data by subtracting
    low-pass components obtained using a Gaussian filter defined by a
    cut-off time period (in days).
    :param tsincr: Array of incremental time series data of shape (ifg.shape, n_epochs).
    :param epochlist: A pyrate.core.shared.EpochList Class instance.
    :param params: Dictionary of PyRate configuration parameters.
    :return: ts_hp: Filtered high frequency time series data; shape (ifg.shape, nepochs).
    """
    log.info('Applying temporal high-pass filter')
    threshold = params[C.TLPF_PTHR]
    cutoff_day = params[C.TLPF_CUTOFF]
    if cutoff_day < 1 or type(cutoff_day) != int:
        raise ValueError(f'tlpf_cutoff must be an integer greater than or '
                         f'equal to 1 day. Value provided = {cutoff_day}')

    # convert cutoff in days to years
    cutoff_yr = cutoff_day / ifc.DAYS_PER_YEAR
    log.info(f'Gaussian temporal filter cutoff is {cutoff_day} days '
             f'({cutoff_yr:.4f} years)')

    intv = np.diff(epochlist.spans)  # time interval for the neighboring epochs
    span = epochlist.spans[: tsincr.shape[2]] + intv / 2  # accumulated time
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


def gaussian_temporal_filter(tsincr: np.ndarray, cutoff: float, span: np.ndarray,
                             thr: int) -> np.ndarray:
    """
    Function to apply a Gaussian temporal low-pass filter to a 1D time-series
    vector for one pixel with irregular temporal sampling.
    :param tsincr: 1D time-series vector to be filtered.
    :param cutoff: filter cutoff in years.
    :param span: 1D vector of cumulative time spans, in years.
    :param thr: threshold for non-NaN values in tsincr.
    :return: ts_lp: Low-pass filtered time series vector.
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


def _kernel(x: np.ndarray, sigma: float) -> np.ndarray:
    """
    Gaussian low-pass filter kernel
    """
    return np.exp(-0.5 * (x / sigma) ** 2)
