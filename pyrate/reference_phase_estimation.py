"""
This is a python implementation of the refphsest.m of pirate.
"""
__author__ = 'Sudipta Basak'
__date_created__ = '22/12/15'

import numpy as np
from pyrate import config as cf


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
    number_ifgs = len(ifgs)
    ref_phs = np.zeros(number_ifgs)
    _validate_ifgs(ifgs)

    # set reference phase as the average of the whole image (recommended)
    if int(params[cf.REF_EST_METHOD]) == 1:
        ifg_phase_data_sum = np.zeros(ifgs[0].shape, dtype=np.float64)
        for i in ifgs:
            ifg_phase_data_sum += i.phase_data
        comp = np.isnan(ifg_phase_data_sum)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
        for n, i in enumerate(ifgs):
            ifgv = np.ravel(i.phase_data, order='F')
            ifgv[comp == 1] = np.nan
            # reference phase
            ref_phs[n] = np.nanmedian(ifgv)
            i.phase_data -= ref_phs[n]

    elif int(params[cf.REF_EST_METHOD]) == 2:
        half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE]/2.0))
        chipsize = 2 * half_chip_size + 1
        thresh = chipsize*chipsize*params[cf.REF_MIN_FRAC]
        for n, i in enumerate(ifgs):
            patch = i.phase_data[
                    refpy - half_chip_size: refpy + half_chip_size + 1,
                    refpx - half_chip_size: refpx + half_chip_size + 1
                    ]
            patch = np.reshape(patch, newshape=(-1, 1), order='F')
            if np.sum(~np.isnan(patch)) < thresh:
                raise ReferencePhaseError('The reference pixel'
                                          'is not in high coherent area!')
            ref_phs[n] = np.nanmedian(patch)
            i.phase_data -= ref_phs[n]
    else:
        raise ReferencePhaseError('No such option. Use refest=1 or 2')

    return ref_phs, ifgs


def _validate_ifgs(ifgs):
    if len(ifgs) < 2:
        raise ReferencePhaseError('Need to provide at least 2 ifgs')


class ReferencePhaseError(Exception):
    """
    Generic class for errors in reference phase estimation.
    """
    pass

if __name__ == "__main__":
    import os
    import shutil
    from subprocess import call

    from pyrate.scripts import run_pyrate
    from pyrate import matlab_mst_kruskal as matlab_mst
    from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR, SYD_TEST_OUT

    # start each full test run cleanly
    shutil.rmtree(SYD_TEST_OUT, ignore_errors=True)

    os.makedirs(SYD_TEST_OUT)

    params = cf.get_config_params(
            os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'pyrate_system_test.conf'))

    call(["python", "pyrate/scripts/run_prepifg.py",
          os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'pyrate_system_test.conf')])

    xlks, ylks, crop = run_pyrate.transform_params(params)

    base_ifg_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

    dest_paths = run_pyrate.get_dest_paths(base_ifg_paths, crop, params, xlks)

    ifg_instance = matlab_mst.IfgListPyRate(datafiles=dest_paths)

    ifgs = ifg_instance.ifgs
    nan_conversion = int(params[cf.NAN_CONVERSION])
    for i in ifgs:
        i.convert_to_mm()
        i.write_modified_phase()
        if nan_conversion and not i.nan_converted:  # nan conversion happens here in networkx mst
            i.convert_to_nans()

    refx, refy = run_pyrate.find_reference_pixel(ifgs, params)

    if params[cf.ORBITAL_FIT] != 0:
        run_pyrate.remove_orbital_error(ifgs, params)

    ref_phs, ifgs = estimate_ref_phase(ifgs, params, refx, refy)
    print ref_phs
