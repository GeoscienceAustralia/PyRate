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
This Python module implements residual orbital corrections for interferograms.
"""
# pylint: disable=invalid-name
import logging
from collections import OrderedDict
from numpy import empty, isnan, reshape, float32, squeeze
from numpy import dot, vstack, zeros, meshgrid
import numpy as np
from numpy.linalg import pinv
# from joblib import Parallel, delayed
from scipy.linalg import lstsq

from pyrate.algorithm import master_slave_ids, get_all_epochs
from pyrate import mst, shared, prepifg
from pyrate.shared import nanmedian, Ifg
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import mpiops

log = logging.getLogger(__name__)

# Orbital correction tasks
#
# TODO: options for multilooking
# 1) do the 2nd stage mlook at prepifg.py/generate up front, then delete in
#    workflow afterward
# 2) refactor prep_ifgs() call to take input filenames and params & generate
#    mlooked versions from that
#    this needs to be more generic to call at any point in the runtime.


# Design notes:
# The orbital correction code is based on MATLAB Pirate, but includes several
# enhancements. Pirate creates sparse arrays for the linear inversion, which
# contain many empty cells. This is unnecessary for the independent method, and
# temporarily wastes potentially a lot of memory.
#
# For the independent method, PyRate makes individual small design matrices and
# corrects the Ifgs one by one. If required in the correction, the offsets
# option adds an extra colum of ones to include in the inversion.
#
# Network method design matrices are mostly empty, and offsets are handled
# differently. Individual design matrices (== independent method DMs) are
# placed in the sparse network design matrix. Offsets are not included in the
# smaller DMs to prevent unwanted cols being inserted. This is why some funcs
# appear to ignore the offset parameter in the networked method. Network DM
# offsets are cols of 1s in a diagonal line on the LHS of the sparse array.

# ORBITAL ERROR correction constants
INDEPENDENT_METHOD = cf.INDEPENDENT_METHOD
NETWORK_METHOD = cf.NETWORK_METHOD

PLANAR = cf.PLANAR
QUADRATIC = cf.QUADRATIC
PART_CUBIC = cf.PART_CUBIC


def orbital_correction(ifgs_or_ifg_paths, params, mlooked=None, offset=True,
                       preread_ifgs=None):
    """
    Removes orbital error from given Ifgs.

    NB: the ifg data is modified in situ, rather than create intermediate files.
    The network method assumes the given ifgs have already been reduced to a
    minimum set from an MST type operation.

    :param ifgs: sequence of Ifg objs to correct
    :param degree: PLANAR, QUADRATIC or PART_CUBIC
    :param method: INDEPENDENT_METHOD or NETWORK_METHOD
    :param mlooked: sequence of multilooked ifgs (must correspond to 'ifgs' arg)
    :param bool offset: True/False to include the constant/offset component
    """
    degree = params[cf.ORBITAL_FIT_DEGREE]
    method = params[cf.ORBITAL_FIT_METHOD]
    # parallel = params[cf.PARALLEL]  # not implemented

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        msg = "Invalid degree of %s for orbital correction" % degree
        raise OrbitalError(msg)

    if method == NETWORK_METHOD:
        if mlooked is None:
            network_correction(ifgs_or_ifg_paths, degree, offset, params,
                               m_ifgs=mlooked,
                               preread_ifgs=preread_ifgs)
        else:
            _validate_mlooked(mlooked, ifgs_or_ifg_paths)
            network_correction(ifgs_or_ifg_paths, degree, offset, params,
                               mlooked, preread_ifgs)

    elif method == INDEPENDENT_METHOD:
        # not running in parallel
        # raises swig object pickle error
        # Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
        #     delayed(_independent_correction)(ifg, degree, offset, params)
        #     for ifg in ifgs)
        for ifg in ifgs_or_ifg_paths:
            independent_correction(ifg, degree, offset, params)
    else:
        msg = "Unknown method: '%s', need INDEPENDENT or NETWORK method"
        raise OrbitalError(msg % method)


def _validate_mlooked(mlooked, ifgs):
    '''
    Basic sanity checking of the multilooked ifgs.
    '''

    if len(mlooked) != len(ifgs):
        msg = "Mismatching # ifgs and # multilooked ifgs"
        raise OrbitalError(msg)

    if not all([hasattr(i, 'phase_data') for i in mlooked]):
        msg = "Mismatching types in multilooked ifgs arg:\n%s" % mlooked
        raise OrbitalError(msg)


def _get_num_params(degree, offset=None):
    '''
    Returns number of model parameters
    '''

    if degree == PLANAR:
        nparams = 2
    elif degree == QUADRATIC:
        nparams = 5
    elif degree == PART_CUBIC:
        nparams = 6
    else:
        msg = "Invalid orbital model degree: %s" % degree
        raise OrbitalError(msg)

    # NB: independent method only, network method handles offsets separately
    if offset is True:
        nparams += 1  # eg. y = mx + offset
    return nparams


def independent_correction(ifg, degree, offset, params):
    """
    Calculates and removes orbital correction from an ifg.

    .. warn:: This will write orbital error corrected phase_data in the ifgs.

    :param ifg: the ifg to remove remove the orbital error from
    :param degree: type of model to use PLANAR, QUADRATIC etc
    :param offset: boolean
    :param params: parameter dictionary
    """
    ifg = shared.Ifg(ifg) if isinstance(ifg, str) else ifg
    if not ifg.is_open:
        ifg.open()
    shared.nan_and_mm_convert(ifg, params)
    # vectorise, keeping NODATA
    vphase = reshape(ifg.phase_data, ifg.num_cells)
    dm = get_design_matrix(ifg, degree, offset)

    # filter NaNs out before getting model
    clean_dm = dm[~isnan(vphase)]
    data = vphase[~isnan(vphase)]
    model = lstsq(clean_dm, data)[0]  # first arg is the model params

    # calculate forward model & morph back to 2D
    if offset:
        fullorb = np.reshape(np.dot(dm[:, :-1], model[:-1]),
                             ifg.phase_data.shape)
    else:
        fullorb = np.reshape(np.dot(dm, model), ifg.phase_data.shape)
    offset_removal = nanmedian(np.ravel(ifg.phase_data - fullorb))
    ifg.phase_data -= (fullorb - offset_removal)
    # set orbfit tags after orbital error correction
    _save_orbital_error_corrected_phase(ifg)
    if ifg.open():
        ifg.close()


def network_correction(ifgs, degree, offset, params, m_ifgs=None,
                       preread_ifgs=None):
    """
    Calculates orbital correction model, removing this from the ifgs.
    .. warn:: This will write orbital error corrected phase_data in the ifgs.

    :param ifgs: interferograms reduced to a minimum tree from prior
        MST calculations
    :param degree: PLANAR, QUADRATIC or PART_CUBIC
    :param offset: True to calculate the model using offsets
    :param params: dict, parameter dictionary
    :param m_ifgs: multilooked orbfit ifgs (sequence must be mlooked
        versions of 'ifgs' arg)
    :param preread_ifgs: dict, parameters dict corresponding to config file
    """
    # pylint: disable=too-many-locals, too-many-arguments
    src_ifgs = ifgs if m_ifgs is None else m_ifgs
    src_ifgs = mst.mst_from_ifgs(src_ifgs)[3]  # use networkx mst

    vphase = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in src_ifgs])
    vphase = squeeze(vphase)

    B = get_network_design_matrix(src_ifgs, degree, offset)

    # filter NaNs out before getting model
    B = B[~isnan(vphase)]
    orbparams = dot(pinv(B, 1e-6), vphase[~isnan(vphase)])

    ncoef = _get_num_params(degree)
    if preread_ifgs:
        temp_ifgs = OrderedDict(sorted(preread_ifgs.items())).values()
        ids = master_slave_ids(get_all_epochs(temp_ifgs))
    else:
        ids = master_slave_ids(get_all_epochs(ifgs))
    coefs = [orbparams[i:i+ncoef] for i in
             range(0, len(set(ids)) * ncoef, ncoef)]

    # create full res DM to expand determined coefficients into full res
    # orbital correction (eg. expand coarser model to full size)

    if preread_ifgs:
        temp_ifg = Ifg(ifgs[0])  # ifgs here are paths
        temp_ifg.open()
        dm = get_design_matrix(temp_ifg, degree, offset=False)
        temp_ifg.close()
    else:
        dm = get_design_matrix(ifgs[0], degree, offset=False)

    # TODO: remove this import from tests suite
    from tests.common import MockIfg
    for i in ifgs:
        # open if not Ifg instance
        if not (isinstance(i, Ifg) or isinstance(i, MockIfg)):  # pragma: no cover
            # are paths
            i = Ifg(i)
            i.open(readonly=False)
            shared.nan_and_mm_convert(i, params)
        _remove_networkx_error(coefs, dm, i, ids, offset)


def _remove_networkx_error(coefs, dm, ifg, ids, offset):
    orb = dm.dot(coefs[ids[ifg.slave]] - coefs[ids[ifg.master]])
    orb = orb.reshape(ifg.shape)
    # offset estimation
    if offset:
        # bring all ifgs to same base level
        orb -= nanmedian(np.ravel(ifg.phase_data - orb))
    ifg.phase_data -= orb  # remove orbital error from the ifg
    # set orbfit tags after orbital error correction
    _save_orbital_error_corrected_phase(ifg)


def _save_orbital_error_corrected_phase(ifg):
    """
    Convenceince function to update metadata and save latest phase after
    orbital fit correction
    Parameters
    ----------
    ifg: Ifg class instance
    """
    # set orbfit tags after orbital error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR, ifc.ORB_REMOVED)
    ifg.write_modified_phase()
    ifg.close()


# TODO: subtract reference pixel coordinate from x and y
def get_design_matrix(ifg, degree, offset, scale=100.0):
    """
    Returns simple design matrix with columns for model parameters.

    :param ifg: interferogram to base the DM on
    :param degree: PLANAR, QUADRATIC or PART_CUBIC
    :param offset: True to include offset cols, otherwise False.
    :param scale: amount to divide cell size by for improving inversion robustness
    """

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        raise OrbitalError("Invalid degree argument")

    # scaling required with higher degree models to help with param estimation
    xsize = ifg.x_size / scale if scale else ifg.x_size
    ysize = ifg.y_size / scale if scale else ifg.y_size

    # mesh needs to start at 1, otherwise first cell resolves to 0 and ignored
    xg, yg = [g+1 for g in meshgrid(range(ifg.ncols), range(ifg.nrows))]
    x = xg.reshape(ifg.num_cells) * xsize
    y = yg.reshape(ifg.num_cells) * ysize

    # TODO: performance test this vs np.concatenate (n by 1 cols)??
    dm = empty((ifg.num_cells, _get_num_params(degree, offset)), dtype=float32)

    # apply positional parameter values, multiply pixel coordinate by cell size
    # to get distance (a coord by itself doesn't tell us distance from origin)
    if degree == PLANAR:
        dm[:, 0] = x
        dm[:, 1] = y
    elif degree == QUADRATIC:
        dm[:, 0] = x**2
        dm[:, 1] = y**2
        dm[:, 2] = x * y
        dm[:, 3] = x
        dm[:, 4] = y
    elif degree == PART_CUBIC:
        dm[:, 0] = x * (y**2)
        dm[:, 1] = x**2
        dm[:, 2] = y**2
        dm[:, 3] = x * y
        dm[:, 4] = x
        dm[:, 5] = y
    if offset is True:
        dm[:, -1] = np.ones(ifg.num_cells)

    return dm


def get_network_design_matrix(ifgs, degree, offset):
    # pylint: disable=too-many-locals
    """
    Returns larger format design matrix for networked error correction.

    The network design matrix includes rows which relate to those of NaN cells.
    :param ifgs: sequence of interferograms
    :param degree: PLANAR, QUADRATIC or PART_CUBIC
    :param offset: True to include offset cols, otherwise False.
    """

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        raise OrbitalError("Invalid degree argument")

    nifgs = len(ifgs)
    if nifgs < 1:
        # can feasibly do correction on a single Ifg/2 epochs
        raise OrbitalError("Invalid number of Ifgs: %s" % nifgs)

    # init sparse network design matrix
    nepochs = len(set(get_all_epochs(ifgs)))

    # no offsets: they are made separately below
    ncoef = _get_num_params(degree)
    shape = [ifgs[0].num_cells * nifgs, ncoef * nepochs]

    if offset:
        shape[1] += nifgs  # add extra block for offset cols

    netdm = zeros(shape, dtype=float32)

    # calc location for individual design matrices
    dates = [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]
    ids = master_slave_ids(dates)
    offset_col = nepochs * ncoef  # base offset for the offset cols
    tmpdm = get_design_matrix(ifgs[0], degree, offset=False)

    # iteratively build up sparse matrix
    for i, ifg in enumerate(ifgs):
        rs = i * ifg.num_cells  # starting row
        m = ids[ifg.master] * ncoef  # start col for master
        s = ids[ifg.slave] * ncoef  # start col for slave
        netdm[rs:rs + ifg.num_cells, m:m + ncoef] = -tmpdm
        netdm[rs:rs + ifg.num_cells, s:s + ncoef] = tmpdm

        # offsets are diagonal cols across the extra array block created above
        if offset:
            netdm[rs:rs + ifg.num_cells, offset_col + i] = 1  # init offset cols

    return netdm


class OrbitalError(Exception):
    """
    Generic class for errors in orbital correction.
    """


def remove_orbital_error(ifgs, params, preread_ifgs=None):
    """
    Wrapper for orbital error removal functionality.

    Parameters
    ----------
    ifgs: list
        list of ifgs or ifg paths
    params: dict
        dict corresponding to config parameters
    preread_ifgs: dict, optional
        dict containing information regarding MPI jobs
    """
    log.info('Calculating orbital error correction')

    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required')
        return

    if preread_ifgs:  # don't check except for mpi tests
        # remove non ifg keys
        _ = [preread_ifgs.pop(k) for k in ['gt', 'epochlist', 'md', 'wkt']]
        # perform some general error/sanity checks
        if mpiops.rank == 0:
            _check_orbital_ifgs(preread_ifgs)

    ifg_paths = [i.data_path for i in ifgs] \
        if isinstance(ifgs[0], Ifg) else ifgs

    mlooked = None

    # mlooking is not necessary for independent correction
    # can use multiple procesing if write_to_disc=True
    if params[cf.ORBITAL_FIT_METHOD] == 2:
        mlooked_dataset = prepifg.prepare_ifgs(
            ifg_paths,
            crop_opt=prepifg.ALREADY_SAME_SIZE,
            xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
            ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
            thresh=params[cf.NO_DATA_AVERAGING_THRESHOLD],
            write_to_disc=False)
        mlooked = [Ifg(m[1]) for m in mlooked_dataset]

        for m in mlooked:
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]
            m.convert_to_nans()
            m.convert_to_mm()

    orbital_correction(ifgs, params, mlooked=mlooked,
                       preread_ifgs=preread_ifgs)


def _check_orbital_ifgs(preread_ifgs):  # pragma: no cover

    ifg_paths = sorted(preread_ifgs.keys())
    # preread_ifgs[i].metadata contains ifg metadata
    flags = [ifc.PYRATE_ORBITAL_ERROR in preread_ifgs[i].metadata
             for i in ifg_paths]
    if all(flags):
        log.info('Skipped orbital correction, ifgs already corrected')
        return
    if (sum(flags) < len(flags)) and (sum(flags) > 0):
        log.debug('Detected mix of corrected and uncorrected '
                  'orbital error in ifgs')

        for i, flag in zip(ifg_paths, flags):
            if flag:
                msg = '{}: prior orbital error correction detected'.format(i)
            else:
                msg = '{}: no orbital correction detected'.format(i)
            log.debug(msg)
            raise OrbitalError(msg)
    log.info('Orbital error status checked')
