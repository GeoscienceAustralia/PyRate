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
This Python module implements a reference phase estimation algorithm
and is based on the function 'refphsest.m' of the Matlab Pirate package.
"""
from __future__ import print_function
import logging
import numpy as np
from joblib import Parallel, delayed
from pyrate import config as cf
from pyrate.shared import nanmedian
from pyrate import ifgconstants as ifc

log = logging.getLogger(__name__)


def estimate_ref_phase(ifgs, params, refpx, refpy):
    """
    xxxxx
    
    :param ifgs: List of interferograms
    :param params: Parameters of the simulation
    :param refpx: Reference pixel X found by ref pixel method
    :param refpy: Reference pixel Y found by ref pixel method
    
    :return ref_phs: Reference phase correction
    :return ifgs: Reference phase data removed from list of interferograms
    """
    _validate_ifgs(ifgs)

    # set reference phase as the average of the whole image (recommended)
    if params[cf.REF_EST_METHOD] == 1:
        ref_phs = est_ref_phase_method1(ifgs, params)

    elif params[cf.REF_EST_METHOD] == 2:
        ref_phs = est_ref_phase_method2(ifgs, params, refpx, refpy)
    else:
        raise ReferencePhaseError('No such option. Use refest=1 or 2')

    for i in ifgs:
        i.meta_data[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
        i.write_modified_phase()
    return ref_phs, ifgs


def est_ref_phase_method2(ifgs, params, refpx, refpy):
    """
    Reference phase estimate method 2. xxxx
    
    :param ifgs: List of interferograms
    :param params: Parameters of the simulation
    :param refpx: Reference pixel X found by ref pixel method
    :param refpy: Reference pixel Y found by ref pixel method
    
    :return ref_phs: Reference phase correction
    :return ifgs: Reference phase data removed from list of interferograms
    """
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    phase_data = [i.phase_data for i in ifgs]
    if params[cf.PARALLEL]:
        ref_phs = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(est_ref_phs_method2)(p, half_chip_size,
                                         refpx, refpy, thresh)
            for p in phase_data)

        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = \
                est_ref_phs_method2(phase_data[n], half_chip_size,
                                    refpx, refpy, thresh)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def est_ref_phs_method2(phase_data, half_chip_size,
                        refpx, refpy, thresh):
    """
    Convenience function for reference phase estimate method 2 parallelisation.
    
    :param phase_data: xxxx
    :param half_chip_size: xxxx
    :param refpx: xxxx
    :param refpy: xxxx
    :param thres: xxxx
    
    :return xxxx
    """
    patch = phase_data[refpy - half_chip_size: refpy + half_chip_size + 1,
                       refpx - half_chip_size: refpx + half_chip_size + 1]
    patch = np.reshape(patch, newshape=(-1, 1), order='F')
    if np.sum(~np.isnan(patch)) < thresh:
        raise ReferencePhaseError('The data window at the reference pixel '
                                  'does not have enough valid observations')
    ref_ph = nanmedian(patch)
    return ref_ph


def est_ref_phase_method1(ifgs, params):
    """
    Reference phase estimate method 1 estimation.

    :param ifgs: List of interferograms or shared.IfgPart class instances
    :param params: Parameter dictionary corresponding to config file

    :return ref_phs: Numpy array of size (nifgs, 1)
    """
    ifg_phase_data_sum = np.zeros(ifgs[0].shape, dtype=np.float64)
    phase_data = [i.phase_data for i in ifgs]
    for ifg in ifgs:
        ifg_phase_data_sum += ifg.phase_data

    comp = np.isnan(ifg_phase_data_sum)  # this is the same as in Matlab
    comp = np.ravel(comp, order='F')  # this is the same as in Matlab
    if params[cf.PARALLEL]:
        log.info("Calculating ref phase using multiprocessing")
        ref_phs = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(est_ref_phs_method1)(p, comp)
            for p in phase_data)
        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        log.info("Calculating ref phase")
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = est_ref_phs_method1(ifg.phase_data, comp)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def est_ref_phs_method1(phase_data, comp):
    """
    Convenience function for ref phs estimate method 1 parallelisation.
    
    :param phase_data: xxxx
    :param comp: xxxx
    
    :return xxxx
    """
    ifgv = np.ravel(phase_data, order='F')
    ifgv[comp == 1] = np.nan
    # reference phase
    ref_ph = nanmedian(ifgv)
    return ref_ph


def _validate_ifgs(ifgs):
    """
    Make sure all interferograms have the same reference phase status.
    """
    if len(ifgs) < 2:
        raise ReferencePhaseError('Need to provide at least 2 ifgs')
    flags = [i.dataset.GetMetadataItem(ifc.REF_PHASE) for i in ifgs]
    if all(flags):
        log.info('Ifgs already reference phase corrected')
        return
    else:
        check_ref_phase_ifgs(ifgs, flags)


def check_ref_phase_ifgs(ifgs, flags):
    """
    Function to check that the reference phase status of all interferograms are the same.
    
    :param ifgs: xxxx
    :param flags: xxxx
    
    :return xxxx    
    """
    count = sum([f == ifc.REF_PHASE_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        log.debug('Detected mix of corrected and uncorrected '
                  'reference phases in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '{}: prior reference phase ' \
                      'correction detected'.format(i.data_path)
            else:
                msg = '{}: no reference phase ' \
                      'correction detected'.format(i.data_path)
            log.debug(msg.format(i.data_path))
            raise ReferencePhaseError(msg)


class ReferencePhaseError(Exception):
    """
    Generic class for errors in reference phase estimation.
    """
    pass
