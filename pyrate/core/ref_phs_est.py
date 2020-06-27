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
# coding: utf-8
"""
This Python module implements a reference phase estimation algorithm.
"""
from joblib import Parallel, delayed
import numpy as np

from pyrate.core import ifgconstants as ifc, config as cf
from pyrate.core.shared import joblib_log_level, nanmedian, Ifg
from pyrate.core import mpiops

from pyrate.core.logger import pyratelogger as log


MASTER_PROCESS = 0


def est_ref_phase_method2(ifg_paths, params, refpx, refpy):
    """
    Reference phase estimation using method 2. Reference phase is the
    median calculated with a patch around the supplied reference pixel.

    :param list ifg_paths: List of interferogram paths or objects.
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

    def _inner(ifg_paths):
        if isinstance(ifg_paths[0], Ifg):
            ifgs = ifg_paths
        else:
            ifgs = [Ifg(ifg_path) for ifg_path in ifg_paths]

        for ifg in ifgs:
            if not ifg.is_open:
                ifg.open(readonly=False)

        phase_data = [i.phase_data for i in ifgs]
        if params[cf.PARALLEL]:
            ref_phs = Parallel(n_jobs=params[cf.PROCESSES],
                               verbose=joblib_log_level(cf.LOG_LEVEL))(
                delayed(_est_ref_phs_method2)(p, half_chip_size, refpx, refpy, thresh)
                for p in phase_data)

            for n, ifg in enumerate(ifgs):
                ifg.phase_data -= ref_phs[n]
        else:
            ref_phs = np.zeros(len(ifgs))
            for n, ifg in enumerate(ifgs):
                ref_phs[n] = _est_ref_phs_method2(phase_data[n], half_chip_size, refpx, refpy, thresh)
                ifg.phase_data -= ref_phs[n]

        for ifg in ifgs:
            _update_phase_metadata(ifg)
            ifg.close()

        return ref_phs
    
    process_ifgs_paths = mpiops.array_split(ifg_paths)
    ref_phs = _inner(process_ifgs_paths)
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


def est_ref_phase_method1(ifg_paths, params):
    """
    Reference phase estimation using method 1. Reference phase is the
    median of the whole interferogram image.

    :param list ifg_paths: List of interferogram paths or objects
    :param dict params: Dictionary of configuration parameters

    :return: ref_phs: Numpy array of reference phase values of size (nifgs, 1)
    :rtype: ndarray
    :return: ifgs: Reference phase data is removed interferograms in place
    """
    def _process_phase_sum(ifg_paths):
        if isinstance(ifg_paths[0], Ifg):
            proc_ifgs = ifg_paths
        else:
            proc_ifgs = [Ifg(ifg_path) for ifg_path in ifg_paths]

        for ifg in proc_ifgs:
            if not ifg.is_open:
                ifg.open(readonly=False)

        ifg_phase_data_sum = np.zeros(proc_ifgs[0].shape, dtype=np.float32)

        for ifg in proc_ifgs:
            ifg_phase_data_sum += ifg.phase_data

        return ifg_phase_data_sum

    def _inner(proc_ifgs, phase_data_sum):
        if isinstance(proc_ifgs[0], Ifg):
            proc_ifgs = proc_ifgs
        else:
            proc_ifgs = [Ifg(ifg_path) for ifg_path in proc_ifgs]

        for ifg in proc_ifgs:
            if not ifg.is_open:
                ifg.open(readonly=False)

        comp = np.isnan(phase_data_sum)
        comp = np.ravel(comp, order='F')

        if params[cf.PARALLEL]:
            log.info("Calculating ref phase using multiprocessing")
            ref_phs = Parallel(n_jobs=params[cf.PROCESSES], verbose=joblib_log_level(cf.LOG_LEVEL))(
                delayed(_est_ref_phs_method1)(p.phase_data, comp) for p in proc_ifgs
            )
            for n, ifg in enumerate(proc_ifgs):
                ifg.phase_data -= ref_phs[n]
        else:
            log.info("Calculating ref phase")
            ref_phs = np.zeros(len(proc_ifgs))
            for n, ifg in enumerate(proc_ifgs):
                ref_phs[n] = _est_ref_phs_method1(ifg.phase_data, comp)
                ifg.phase_data -= ref_phs[n]

        for ifg in proc_ifgs:
            _update_phase_metadata(ifg)
            ifg.close()

        return ref_phs

    process_ifg_paths = mpiops.array_split(ifg_paths)
    ifg_phase_data_sum = mpiops.comm.allreduce(_process_phase_sum(process_ifg_paths), mpiops.sum0_op)
    ref_phs = _inner(process_ifg_paths, ifg_phase_data_sum)

    return ref_phs


def _est_ref_phs_method1(phase_data, comp):
    """
    Convenience function for ref phs estimate method 1 parallelisation
    """
    ifgv = np.ravel(phase_data, order='F')
    ifgv[comp == 1] = np.nan
    return nanmedian(ifgv)


def _update_phase_metadata(ifg):
    ifg.meta_data[ifc.PYRATE_REF_PHASE] = ifc.REF_PHASE_REMOVED
    ifg.write_modified_phase()
    log.debug(f"Reference phase corrected for {ifg.data_path}")


class ReferencePhaseError(Exception):
    """
    Generic class for errors in reference phase estimation.
    """
    pass
