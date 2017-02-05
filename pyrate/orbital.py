'''
Module for calculating orbital correction for interferograms.

Created on 31/3/13

.. codeauthor:: Ben Davies, Sudipta Basak
'''

from numpy import empty, isnan, reshape, float32, squeeze
from numpy import dot, vstack, zeros, meshgrid
# from joblib import Parallel, delayed
import numpy as np
from scipy.linalg import lstsq
from numpy.linalg import pinv

from pyrate.algorithm import master_slave_ids, get_all_epochs, get_epoch_count
from pyrate import mst, shared
from pyrate.shared import nanmedian
from pyrate import config as cf
from pyrate import ifgconstants as ifc


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


def orbital_correction(ifgs_or_ifg_paths, params, mlooked=None, offset=True):
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
    parallel = params[cf.PARALLEL]  # not implemented

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        msg = "Invalid degree of %s for orbital correction" % degree
        raise OrbitalError(msg)

    if method == NETWORK_METHOD:
        ifgs = ifgs_or_ifg_paths
        if mlooked is None:
            _network_correction(ifgs, degree, offset)
        else:
            _validate_mlooked(mlooked, ifgs)
            _network_correction(ifgs, degree, offset, mlooked)

    elif method == INDEPENDENT_METHOD:
        # not running in parallel
        # raises swig object pickle error
        # Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
        #     delayed(_independent_correction)(ifg, degree, offset, params)
        #     for ifg in ifgs)
        for ifg in ifgs_or_ifg_paths:
            _independent_correction(ifg, degree, offset, params)
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


def get_num_params(degree, offset=None):
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


def _independent_correction(ifg, degree, offset, params):
    """
    Calculates and removes orbital correction from an ifg.

    .. warn:: Changes are made in place to the Ifg obj.

    :param ifg: the ifg to remove remove the orbital error from
    :param degree: type of model to use PLANAR, QUADRATIC etc
    :param offset: boolean
    :param params: parameter dictionary
    """
    ifg = shared.Ifg(ifg) if isinstance(ifg, str) else ifg
    if not ifg.is_open:
        ifg.open()
    shared.nan_and_mm_convert(ifg, params)
    vphase = reshape(ifg.phase_data, ifg.num_cells)  # vectorise, keeping NODATA
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
    save_orbital_error_corrected_phase(ifg)
    if ifg.open():
        ifg.close()


def _network_correction(ifgs, degree, offset, m_ifgs=None):
    """
    Calculates orbital correction model, removing this from the ifgs.
    .. warn:: This does in-situ modification of phase_data in the ifgs.

    :param ifgs: interferograms reduced to a minimum tree from prior MST calculations
    :param degree: PLANAR, QUADRATIC or PART_CUBIC
    :param offset: True to calculate the model using offsets
    :param m_ifgs: multilooked ifgs (sequence must be mlooked versions of 'ifgs' arg)
    """
    # get DM & filter out NaNs
    src_ifgs = ifgs if m_ifgs is None else m_ifgs
    src_ifgs = mst.mst_from_ifgs(src_ifgs)[3]  # use networkx mst

    vphase = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in src_ifgs])
    vphase = squeeze(vphase)

    B = get_network_design_matrix(src_ifgs, degree, offset)

    # filter NaNs out before getting model
    B = B[~isnan(vphase)]
    obsv = vphase[~isnan(vphase)]
    orbparams = dot(pinv(B, 1e-6), obsv)

    ncoef = get_num_params(degree)
    ids = master_slave_ids(get_all_epochs(ifgs))
    coefs = [orbparams[i:i+ncoef] for i in
             range(0, len(set(ids)) * ncoef, ncoef)]

    # create full res DM to expand determined coefficients into full res orbital
    # correction (eg. expand coarser model to full size)
    ifg_shape = ifgs[0].shape
    dm = get_design_matrix(ifgs[0], degree, offset=False)

    for i in ifgs:
        orb = dm.dot(coefs[ids[i.slave]] - coefs[ids[i.master]])
        orb = orb.reshape(ifg_shape)

        # offset estimation
        if offset:
            tmp = np.ravel(i.phase_data - orb)
            # bring all ifgs to same base level
            orb -= nanmedian(tmp)

        i.phase_data -= orb  # remove orbital error from the ifg
        # set orbfit tags after orbital error correction
        save_orbital_error_corrected_phase(i)


def save_orbital_error_corrected_phase(ifg):
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
    dm = empty((ifg.num_cells, get_num_params(degree, offset)), dtype=float32)

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
    nepochs = get_epoch_count(ifgs)
    ncoef = get_num_params(degree)  # no offsets: they are made separately below
    shape = [ifgs[0].num_cells * nifgs, ncoef * nepochs]

    if offset:
        shape[1] += nifgs # add extra block for offset cols

    netdm = zeros(shape, dtype=float32)

    # calc location for individual design matrices
    dates = [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]
    ids = master_slave_ids(dates)
    offset_col = nepochs * ncoef # base offset for the offset cols
    tmpdm = get_design_matrix(ifgs[0], degree, offset=False)

    # iteratively build up sparse matrix
    for i, ifg in enumerate(ifgs):
        rs = i * ifg.num_cells # starting row
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

    pass
