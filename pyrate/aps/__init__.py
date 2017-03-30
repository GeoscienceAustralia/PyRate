import logging
import os
from copy import deepcopy
from collections import OrderedDict
import numpy as np

from pyrate import config as cf, mpiops, shared
from pyrate.algorithm import get_epochs
from pyrate.aps.spatial import spatial_low_pass_filter
from pyrate.aps.temporal import tlpfilter
from pyrate.scripts.postprocessing import assemble_tiles
from pyrate.shared import Ifg
from pyrate import ifgconstants as ifc
from pyrate.timeseries import time_series

log = logging.getLogger(__name__)


def wrap_spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs):
    """
    Just a wrapper for spatio-tempotal-filter so it can be tested.
    See docstring for spatio_temporal_filter.
    Required due to matlab and python using different mst implementations
    """
    if not params[cf.APSEST]:
        log.info('APS correction not required.')
        return
    tsincr = calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles)

    ifg = Ifg(ifg_paths[0])  # just grab any for parameters in slpfilter
    ifg.open()
    spatio_temporal_filter(tsincr, ifg, params, preread_ifgs)
    ifg.close()


def spatio_temporal_filter(tsincr, ifg, params, preread_ifgs):
    """
    Applies spatio-temportal (aps) filter and save the interferograms.
    
    A first step is to compute the time series using the non-smooth SVD method.
    This is followed by temporal and spatial filters respectively.
     
    Parameters
    ----------
    tsincr: ndarray
        time series array of size (ifg.shape, nvels)
    ifg: shated.Ifg instance
        list of ifg paths
    params: dict
        parameters dict
    tiles: list
        list of shared.Tile class instances
    preread_ifgs: dict
        dictionary of {ifgpath:shared.PrereadIfg class instances}
    """
    epochlist = mpiops.run_once(get_epochs, preread_ifgs)[0]
    tsfilt_incr = mpiops.run_once(tlpfilter, tsincr, epochlist, params)
    ts_lp = tsincr - tsfilt_incr
    ts_aps = mpiops.run_once(spatial_low_pass_filter, ts_lp, ifg, params)
    tsincr -= ts_aps
    mpiops.run_once(ts_to_ifgs, tsincr, preread_ifgs)
    

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

    nvels = None
    for t in process_tiles:
        log.info('Calculating time series for tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir,
                                        'mst_mat_{}.npy'.format(t.index)))
        # don't save this time_series calc, this is done as part of aps filter
        tsincr = time_series(ifg_parts, new_params, vcmt=None, mst=mst_tile)[0]
        np.save(file=os.path.join(output_dir, 'tsincr_aps_{}.npy'.format(
            t.index)), arr=tsincr)
        nvels = tsincr.shape[2]

    nvels = mpiops.comm.bcast(nvels, root=0)
    # need to assemble tsincr from all processes
    tsincr_g = mpiops.run_once(_assemble_tsincr, ifg_paths, params,
                               preread_ifgs, tiles, nvels)
    log.info('Finished calculating time series for spatio-temporal filter')
    return tsincr_g


def _assemble_tsincr(ifg_paths, params, preread_ifgs, tiles, nvels):
    shape = preread_ifgs[ifg_paths[0]].shape + (nvels,)
    tsincr_g = np.empty(shape=shape, dtype=np.float32)
    for i in range(nvels):
        for n, t in enumerate(tiles):
            assemble_tiles(i, n, t, tsincr_g[:, :, i], params[cf.TMPDIR],
                           'tsincr_aps')
    return tsincr_g


def ts_to_ifgs(ts, preread_ifgs):
    """
    Function that takes a time series and converts into ifg phase data
    Parameters
    ----------
    ts: ndarray
        time series nd array (ifg.shape, nvels)
    preread_ifgs: dict
        dict with ifg basic informations
    
    Saves ifgs on disc after conversion.
    """
    log.info('Converting time series to ifgs')
    ifgs = list(OrderedDict(sorted(preread_ifgs.items())).values())
    epochlist, n = get_epochs(ifgs)
    index_master, index_slave = n[:len(ifgs)], n[len(ifgs):]
    for i in range(len(ifgs)):
        phase = np.sum(ts[:, :, index_master[i]: index_slave[i]-1], axis=2)
        _save_aps_corrected_phase(ifgs[i].path, phase)


def _save_aps_corrected_phase(ifg_path, phase):
    """
    Update metadata and save latest phase after
    spatio-temporal filter (aps) correction
    Parameters
    ----------
    ifg_path: str
        path to ifg
    phase: ndarray
        ndarray corresponding to ifg phase data of shape ifg.shape
    """

    ifg = Ifg(ifg_path)
    ifg.open(readonly=False)
    ifg.phase_data[~np.isnan(ifg.phase_data)] = \
        phase[~np.isnan(ifg.phase_data)]
    # set aps tags after aps error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_APS_ERROR, ifc.APS_REMOVED)
    ifg.write_modified_phase()
    ifg.close()
