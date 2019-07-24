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
# coding: utf-8
"""
This Python module implements a reference phase estimation algorithm
and is based on the function 'refphsest.m' of the Matlab Pirate package.
"""
import logging
import numpy as np
from joblib import Parallel, delayed
from pyrate import config as cf
from pyrate.shared import nanmedian
from pyrate import ifgconstants as ifc

log = logging.getLogger(__name__)


def estimate_ref_phase(ifgs, params, refpx, refpy):
    """
    Wrapper function for reference phase estimation and interferogram
    correction.

    :param list ifgs: List of interferogram objects
    :param dict params: Dictionary of configuration parameters
    :param int refpx: Reference pixel X found by ref pixel method
    :param int refpy: Reference pixel Y found by ref pixel method

    :return: ref_phs: Numpy array of reference phase values of size (nifgs, 1)
    :rtype: ndarray
    :return: ifgs: Reference phase data is removed interferograms in place
    """
    _check_ref_phs_ifgs(ifgs)
    # set reference phase as the average of the whole image (recommended)
    if params[cf.REF_EST_METHOD] == 1:
        ref_phs = est_ref_phase_method1(ifgs, params)

    elif params[cf.REF_EST_METHOD] == 2:
        ref_phs = est_ref_phase_method2(ifgs, params, refpx, refpy)
    else:
        raise ReferencePhaseError('No such option. Use refest=1 or 2')

    for i in ifgs:
        i.meta_data[ifc.PYRATE_REF_PHASE] = ifc.REF_PHASE_REMOVED
        i.write_modified_phase()
    return ref_phs, ifgs


def est_ref_phase_method2(ifgs, params, refpx, refpy):
    """
    Reference phase estimation using method 2. Reference phase is the
    median calculated with a patch around the supplied reference pixel.

    :param list ifgs: List of interferogram objects
    :param dict params: Dictionary of configuration parameters
    :param int refpx: Reference pixel X found by ref pixel method
    :param int refpy: Reference pixel Y found by ref pixel method

    :return: ref_phs: Numpy array of reference phase values of size (nifgs, 1)
    :rtype: ndarray
    :return: ifgs: Reference phase data is removed interferograms in place
    """
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    phase_data = [i.phase_data for i in ifgs]
    if params[cf.PARALLEL]:
        ref_phs = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(_est_ref_phs_method2)(p, half_chip_size,
                                          refpx, refpy, thresh)
            for p in phase_data)

        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = \
                _est_ref_phs_method2(phase_data[n], half_chip_size,
                                     refpx, refpy, thresh)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def _est_ref_phs_method2(phase_data, half_chip_size, refpx, refpy, thresh):
    """
    Convenience function for ref phs estimate method 2 parallelisation
    """
    patch = phase_data[refpy - half_chip_size: refpy + half_chip_size + 1,
                       refpx - half_chip_size: refpx + half_chip_size + 1]
    patch = np.reshape(patch, newshape=(-1, 1), order='F')
    nanfrac = np.sum(~np.isnan(patch))
    if nanfrac < thresh:
        raise ReferencePhaseError('The data window at the reference pixel '
                                  'does not have enough valid observations. '
                                  'Actual = {}, Threshold = {}.'.format(
                                          nanfrac, thresh))
    ref_ph = nanmedian(patch)
    return ref_ph


def est_ref_phase_method1(ifgs, params):
    """
    Reference phase estimation using method 1. Reference phase is the
    median of the whole interferogram image.

    :param list ifgs: List of interferogram objects
    :param dict params: Dictionary of configuration parameters

    :return: ref_phs: Numpy array of reference phase values of size (nifgs, 1)
    :rtype: ndarray
    :return: ifgs: Reference phase data is removed interferograms in place
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
            delayed(_est_ref_phs_method1)(p, comp)
            for p in phase_data)
        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        log.info("Calculating ref phase")
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = _est_ref_phs_method1(ifg.phase_data, comp)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def _est_ref_phs_method1(phase_data, comp):
    """
    Convenience function for ref phs estimate method 1 parallelisation
    """
    ifgv = np.ravel(phase_data, order='F')
    ifgv[comp == 1] = np.nan
    return nanmedian(ifgv)


def _check_ref_phs_ifgs(ifgs):
    """
    Convenience function to check the ref phase status of all ifgs
    """
    if len(ifgs) < 2:
        raise ReferencePhaseError('Need to provide at least 2 ifgs')
    # The following code is duplicated from shared.check_correction_status
    flags = [True if i.dataset.GetMetadataItem(ifc.PYRATE_REF_PHASE)
             else False for i in ifgs]

    if sum(flags) == len(flags):
        log.info('Skipped reference phase estimation, ifgs already corrected')
        return True
    elif (sum(flags) < len(flags)) and (sum(flags) > 0):
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
    else:  # count == 0
        log.info('Estimating and removing reference phase')
        return False


class ReferencePhaseError(Exception):
    """
    Generic class for errors in reference phase estimation.
    """
    pass
