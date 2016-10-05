"""
This is a python implementation of the refphsest.m of pirate.
"""
__author__ = 'Sudipta Basak'
__date_created__ = '22/12/15'

import numpy as np
import parmap
from pyrate import config as cf
from pyrate.shared import nanmedian, write_msg
from pyrate import ifgconstants as ifc
import logging


def estimate_ref_phase(ifgs, params, refpx, refpy):
    """
    :param ifgs: list of interferrograms
    :param params: parameters of the simulation
    :param refpx: reference pixel found by ref pixel method
    :param refpy: reference pixel found by ref pixel method
    :returns:
        :ref_phs: reference phase correction
        :ifgs: reference phase data removed list of ifgs
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
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    phase_data = [i.phase_data for i in ifgs]
    if params[cf.PARALLEL]:
        ref_phs = parmap.map(est_ref_phase_method2_multi, phase_data,
                                 half_chip_size, refpx, refpy, thresh)
        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = \
                est_ref_phase_method2_multi(phase_data[n], half_chip_size,
                                            refpx, refpy, thresh)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def est_ref_phase_method2_multi(phase_data, half_chip_size,
                                refpx, refpy, thresh):
    patch = phase_data[
            refpy - half_chip_size: refpy + half_chip_size + 1,
            refpx - half_chip_size: refpx + half_chip_size + 1
            ]
    patch = np.reshape(patch, newshape=(-1, 1), order='F')
    if np.sum(~np.isnan(patch)) < thresh:
        raise ReferencePhaseError('The reference pixel'
                                  'is not in high coherent area!')
    ref_ph = nanmedian(patch)
    return ref_ph


def est_ref_phase_method1(ifgs, params):
    ifg_phase_data_sum = np.zeros(ifgs[0].shape, dtype=np.float64)
    # TODO: revisit as this will likely hit memory limit in NCI
    phase_data = [i.phase_data for i in ifgs]
    for ifg in ifgs:
        ifg_phase_data_sum += ifg.phase_data

    comp = np.isnan(ifg_phase_data_sum)  # this is the same as in Matlab
    comp = np.ravel(comp, order='F')  # this is the same as in Matlab
    if params[cf.PARALLEL]:
        print "ref phase calculation using multiprocessing"
        ref_phs = parmap.map(est_ref_phase_method1_multi, phase_data, comp)
        for n, ifg in enumerate(ifgs):
            ifg.phase_data -= ref_phs[n]
    else:
        print "ref phase calculation in serial"
        ref_phs = np.zeros(len(ifgs))
        for n, ifg in enumerate(ifgs):
            ref_phs[n] = est_ref_phase_method1_multi(ifg.phase_data, comp)
            ifg.phase_data -= ref_phs[n]
    return ref_phs


def est_ref_phase_method1_multi(phase_data, comp):
    ifgv = np.ravel(phase_data, order='F')
    ifgv[comp == 1] = np.nan
    # reference phase
    ref_ph = nanmedian(ifgv)
    return ref_ph


def _validate_ifgs(ifgs):
    if len(ifgs) < 2:
        raise ReferencePhaseError('Need to provide at least 2 ifgs')
    flags = [i.dataset.GetMetadataItem(ifc.REF_PHASE) for i in ifgs]
    if all(flags):
        write_msg('Ifgs already reference phase corrected')
        return
    else:
        check_ref_phase_ifgs(ifgs, flags)


def check_ref_phase_ifgs(ifgs, flags):
    count = sum([f == ifc.REF_PHASE_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected '
                      'reference phases in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '%s: prior reference phase correction detected'
            else:
                msg = '%s: no reference phase correction detected'
            logging.debug(msg % i.data_path)

        raise ReferencePhaseError(msg)


class ReferencePhaseError(Exception):
    """
    Generic class for errors in reference phase estimation.
    """
    pass
