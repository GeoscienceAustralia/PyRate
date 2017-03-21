import os
import logging
from copy import deepcopy
import numpy as np
from numpy import isnan
from pyrate.timeseries import time_series
from pyrate import config as cf
from pyrate import mpiops
from pyrate import shared
from pyrate.algorithm import get_epochs
from pyrate.scripts.postprocessing import assemble_tiles

log = logging.getLogger(__name__)


def spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs):
    log.info('Calculating time series for spatio-temporal filter application')
    if not params[cf.APSEST]:
        log.info('APS correction not required.')
        return
    temporal_low_pass_filter(ifg_paths, params, tiles, preread_ifgs)


def temporal_low_pass_filter(ifg_paths, params, tiles, preread_ifgs):
    log.info('Applying the temporal filter')
    # copy params temporarily
    new_params = deepcopy(params)
    new_params[cf.TIME_SERIES_METHOD] = 2  # use SVD method

    epochlist = mpiops.run_once(get_epochs, preread_ifgs)[0]
    process_tiles = mpiops.array_split(tiles)

    output_dir = params[cf.TMPDIR]
    for t in process_tiles:
        log.info('Calculating time series for tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir,
                                        'mst_mat_{}.npy'.format(t.index)))
        # don't save this time_series calc, this is done as part of aps filter
        tsincr = time_series(ifg_parts, new_params, vcmt=None, mst=mst_tile)[0]
        np.save(file=os.path.join(output_dir, 'tsincr_aps_{}.npy'.format(
            t.index)), arr=tsincr)

    # need to assemble tsincr from all processes
    shape = preread_ifgs[ifg_paths[0]].shape + (tsincr.shape[2], )
    tsincr_g = np.empty(shape=shape, dtype=np.float32)

    for i in range(tsincr.shape[2]):
        for n, t in enumerate(tiles):
            assemble_tiles(i, n, t, tsincr_g[:, :, i], params[cf.TMPDIR],
                           'tsincr_aps')

    ts_hp = tsincr_g - tlpfilter(tsincr_g, epochlist, params)
    log.info('Finished applying temporal low pass filter')


def tlpfilter(tsincr, epochlist, params):
    """temporal low pass filter"""
    nanmat = ~isnan(tsincr)
    tsfilt_incr = np.empty_like(tsincr, dtype=np.float32) * np.nan
    intv = np.diff(epochlist.spans)  # time interval for the neighboring epoch
    span = epochlist.spans[: tsincr.shape[2]] + intv/2  # accumulated time
    rows, cols = tsincr.shape[:2]

    for i in range(rows):
        for j in range(cols):
            sel = np.nonzero(nanmat[i, j, :])[0]  # don't select if nan
            m = len(sel)
            if m >= params[cf.TLPF_PTHR]:
                for k in range(m):
                    yr = span[sel] - span[sel[k]]
                    if params[cf.TLPF_METHOD] == 1:  # gaussian filter
                        wgt = np.exp(-(yr/params[cf.TLPF_CUTOFF]) ** 2/2)
                    elif params[cf.TLPF_METHOD] == 2:  # triangular filter
                        wgt = params[cf.TLPF_CUTOFF] - abs(yr)
                        wgt[wgt < 0] = 0
                    else:  # mean filter
                        wgt = np.ones(m)
                    wgt /= np.sum(wgt)
                    tsfilt_incr[i, j, sel[k]] = np.sum(tsincr[i, j, sel] * wgt)
    return tsfilt_incr
