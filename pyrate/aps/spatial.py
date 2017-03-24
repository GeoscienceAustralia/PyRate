import numpy as np
from pyrate import config as cf
from pyrate.vcm import cvd_from_phase


def slpfilter(ts_hp, params, shape):

    if np.sum(np.isnan(ts_hp)) == 0:  # if there is no nans in the data
        return ts_hp

    if params[cf.SLPF_CUTOFF] == 0:
        maxvar, alpha = cvd_from_phase(ts_hp, ifg, calc_alpha=True)

