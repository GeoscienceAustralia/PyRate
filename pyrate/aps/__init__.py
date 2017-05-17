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
This Python module implements a spatio-temporal filter method
for correcting interferograms for atmospheric phase screen (APS)
signals.
"""

import logging
import os
from copy import deepcopy
from collections import OrderedDict
import numpy as np

from pyrate import config as cf, mpiops, shared
from pyrate.algorithm import get_epochs
from pyrate.aps.spatial import spatial_low_pass_filter
from pyrate.aps.temporal import temporal_low_pass_filter
from pyrate.scripts.postprocessing import _assemble_tiles
from pyrate.shared import Ifg
from pyrate import ifgconstants as ifc
from pyrate.timeseries import time_series

log = logging.getLogger(__name__)


def _wrap_spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs):
    """
    A wrapper for the spatio-temporal filter so it can be tested.
    See docstring for spatio_temporal_filter.
    Required due to differences between Matlab and Python MST
    implementations.
    """
    if not params[cf.APSEST]:
        log.info('APS correction not required.')
        return

    # perform some checks on existing ifgs
    log.info('Checking APS correction status')
    if mpiops.run_once(shared.check_correction_status, preread_ifgs,
                       ifc.PYRATE_APS_ERROR):
        log.info('Finished APS correction')
        return  # return if True condition returned

    tsincr = _calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles)

    ifg = Ifg(ifg_paths[0])  # just grab any for parameters in slpfilter
    ifg.open()
    spatio_temporal_filter(tsincr, ifg, params, preread_ifgs)
    ifg.close()


def spatio_temporal_filter(tsincr, ifg, params, preread_ifgs):
    """
    Applies a spatio-temporal filter to remove the atmospheric phase screen
    (APS) and saves the corrected interferograms. Before performing this step,
    the time series must be computed using the SVD method. This function then
    performs temporal and spatial filtering.

    :param ndarray tsincr: incremental time series array of size
                (ifg.shape, nepochs-1)
    :param list ifg: List of pyrate.shared.Ifg class objects.
    :param dict params: Dictionary of configuration parameter
    :param list tiles: List of pyrate.shared.Tile class objects
    :param dict preread_ifgs: Dictionary of shared.PrereadIfg class instances

    :return: None, corrected interferograms are saved to disk
    """
    epochlist = mpiops.run_once(get_epochs, preread_ifgs)[0]
    ts_lp = mpiops.run_once(temporal_low_pass_filter, tsincr, epochlist, params)
    ts_hp = tsincr - ts_lp
    ts_aps = mpiops.run_once(spatial_low_pass_filter, ts_hp, ifg, params)
    tsincr -= ts_aps

    mpiops.run_once(_ts_to_ifgs, tsincr, preread_ifgs)


def _calc_svd_time_series(ifg_paths, params, preread_ifgs, tiles):
    """
    Helper function to obtain time series for spatio-temporal filter
    using SVD method
    """
    # Is there other existing functions that can perform this same job?
    log.info('Calculating time series via SVD method for '
             'spatio-temporal filter')
    # copy params temporarily
    new_params = deepcopy(params)
    new_params[cf.TIME_SERIES_METHOD] = 2  # use SVD method

    process_tiles = mpiops.array_split(tiles)
    output_dir = params[cf.TMPDIR]

    nvels = None
    for t in process_tiles:
        log.info('Calculating time series for tile {} during aps '
                 'correction'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir,
                                        'mst_mat_{}.npy'.format(t.index)))
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
    """
    Helper function to reconstruct time series images from tiles
    """
    shape = preread_ifgs[ifg_paths[0]].shape + (nvels,)
    tsincr_g = np.empty(shape=shape, dtype=np.float32)
    for i in range(nvels):
        for n, t in enumerate(tiles):
            _assemble_tiles(i, n, t, tsincr_g[:, :, i], params[cf.TMPDIR],
                            'tsincr_aps')
    return tsincr_g


def _ts_to_ifgs(tsincr, preread_ifgs):
    """
    Function that converts an incremental displacement time series into
    interferometric phase observations. Used to re-construct an interferogram
    network from a time series.

    :param ndarray tsincr: incremental time series array of size
                (ifg.shape, nepochs-1)
    :param dict preread_ifgs: Dictionary of shared.PrereadIfg class instances

    :return: None, interferograms are saved to disk
    """
    log.info('Converting time series to ifgs')
    ifgs = list(OrderedDict(sorted(preread_ifgs.items())).values())
    _, n = get_epochs(ifgs)
    index_master, index_slave = n[:len(ifgs)], n[len(ifgs):]
    for i, ifg in enumerate(ifgs):
        phase = np.sum(tsincr[:, :, index_master[i]: index_slave[i]], axis=2)
        _save_aps_corrected_phase(ifg.path, phase)


def _save_aps_corrected_phase(ifg_path, phase):
    """
    Save (update) interferogram metadata and phase data after
    spatio-temporal filter (APS) correction.
    """
    ifg = Ifg(ifg_path)
    ifg.open(readonly=False)
    ifg.phase_data[~np.isnan(ifg.phase_data)] = \
        phase[~np.isnan(ifg.phase_data)]
    # set aps tags after aps error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_APS_ERROR, ifc.APS_REMOVED)
    ifg.write_modified_phase()
    ifg.close()
