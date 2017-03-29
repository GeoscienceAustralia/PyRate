import logging
import os
from copy import deepcopy

import numpy as np

from pyrate import config as cf, mpiops, shared
from pyrate.algorithm import get_epochs
from pyrate.aps.spatial import spatial_low_pass_filter
from pyrate.aps.temporal import tlpfilter
from pyrate.scripts.postprocessing import assemble_tiles
from pyrate.shared import Ifg
from pyrate.timeseries import time_series

log = logging.getLogger(__name__)


def spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs):
    if not params[cf.APSEST]:
        log.info('APS correction not required.')
        return

    tsincr = calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles)
    
    epochlist = mpiops.run_once(get_epochs, preread_ifgs)[0]
    tsfilt_incr = mpiops.run_once(tlpfilter, tsincr, epochlist, params)

    ts_lp = tsincr - tsfilt_incr
    ifg = Ifg(ifg_paths[0])  # just grab any for parameters in slpfilter
    ifg.open()

    ts_aps = spatial_low_pass_filter(ts_lp, ifg, params)
    ifg.close()
    tsincr -= ts_aps

    # need ts2ifgs equivalent here
    # save ifgs
    

def calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles):
    """
    Time series inversion with no smoothing (svd method)
    This is the matlab tsinvnosm.m equivalent.
    
    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict
    preread_ifgs: dict
        prepread ifgs dict
    tiles: list
        list of shared.Tile class instances
    
    Returns
    -------
    tsincr_g: ndarray
        non smooth (svd method) time series of shape (ifg.shape, nvels)
    
    """
    log.info('Calculating time series without smoothing for '
             'spatio-temporal filter')
    # copy params temporarily
    new_params = deepcopy(params)
    new_params[cf.TIME_SERIES_METHOD] = 2  # use SVD method

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
    tsincr_g = mpiops.run_once(_assemble_tsincr, ifg_paths, params,
                               preread_ifgs, tiles, tsincr.shape[2])
    log.info('Finished calculating time series for spatio-temporal filter')
    return tsincr_g


def _assemble_tsincr(ifg_paths, params, preread_ifgs, tiles, nvelpar):
    shape = preread_ifgs[ifg_paths[0]].shape + (nvelpar,)
    tsincr_g = np.empty(shape=shape, dtype=np.float32)
    for i in range(nvelpar):
        for n, t in enumerate(tiles):
            assemble_tiles(i, n, t, tsincr_g[:, :, i], params[cf.TMPDIR],
                           'tsincr_aps')
    return tsincr_g


def ts_to_ifgs(ts, params):
    pass
