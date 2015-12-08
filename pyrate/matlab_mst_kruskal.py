__author__ = 'Sudipta Basak'
__date_created__ = '7/12/15'
import numpy as np

from pyrate.tests.common import sydney_data_setup_ifg_file_list
from pyrate.tests.common import sydney_data_setup
from algorithm import get_epochs_and_n


class IfgList(object):
    """
    copy of matlab ifglist in getnml.m.
    1. Please note that we don't need BaseT unless we are using variance in ifg
    data as cost.
    2.
    """

    def __init__(self, datafiles=None):
        self.nml = sydney_data_setup_ifg_file_list(datafiles)
        self.ifgs = sydney_data_setup(datafiles)
        self.id = range(len(self.nml))
        self.base_T = None
        self.max_var = np.zeros_like(self.id)
        self.alpha = np.zeros_like(self.id)
        self.nan_frac = np.zeros_like(self.id)
        self.master_num = None
        self.slave_num = None
        self.n = None

    def reshape_n(self, n):
        n_reshaped = np.reshape(n, newshape=(2, len(self.id)))
        self.master_num = n_reshaped[0, :]
        self.slave_num = n_reshaped[1, :]
        self.n = n

    def update_nan_frac(self):
        self.nan_frac = [i.nan_fraction for i in self.ifgs]


def get_nml(datafiles=None, nan_conversion=False, prefix_len=4):
    """
    A reproduction of getnml.m, the matlab function in pi-rate.
    Note: the matlab version tested does not have nan's.
    """
    ifg_list_ = IfgList(datafiles)
    epoch_list_, n = get_epochs_and_n(ifg_list_.ifgs)
    ifg_list_.reshape_n(n)
    if nan_conversion:
        ifg_list_.update_nan_frac()  # turn on for nan conversion
    return ifg_list_, epoch_list_


def sort_list(ifg_list_):
    dtype = [('id', int), ('master', int), ('slave', int), ('nan_frac', float)]
    sort_list = map(lambda i, m, s, n: (i, m, s, n),
                    ifg_list_.id, ifg_list_.master_num,
                    ifg_list_.slave_num, ifg_list_.nan_frac)

    np.set_printoptions(precision=3, suppress=True)
    sort_list = np.array(sort_list, dtype=dtype)

    return np.sort(sort_list, order=['nan_frac'])


def matlab_mst_kruskal(ifg_list):
    """
    This is an implementation of the pi-rate mst_kruskal.m
    :param ifg_list:
    :return:
    """

    dtype = [('id', int), ('master', int), ('slave', int), ('nan_frac', float)]
    no_ifgs = len(ifg_list.master_num)
    no_images = max(max(ifg_list.master_num), max(ifg_list.slave_num))
    ifg_sorted = sort_list(ifg_list)

    # add one to ensure index number + 1
    connect = np.eye(no_images + 1)

    mst_list = []

    for i in range(no_ifgs):
        master = ifg_sorted[i][1]
        slave = ifg_sorted[i][2]
        loc_master = np.where(connect[:, master] == 1)[0][0]
        loc_slave = np.where(connect[:, slave] == 1)[0][0]

        if loc_master != loc_slave:
            mst_list.append(ifg_sorted[i])
            connect[loc_master, :] = connect[loc_master, :] + \
                                     connect[loc_slave, :]
            connect = np.delete(connect, loc_slave, axis=0)

    mst_list = np.array(mst_list, dtype=dtype)
    mst_list = np.sort(mst_list, order=['id'])

    return [i[0] for i in mst_list]


def matlab_mst(ifg_list, ifg_list_mst_id, p_thresh_hold=1):
    """
    This is an implementation of matlab/pirate make_mstmat.m.
    """
    data_stack = np.array([i.phase_data for i in ifg_list.ifgs],
                          dtype=np.float32)
    nan_ifg = np.isnan(data_stack)
    mst_mat = np.zeros_like(nan_ifg, dtype=np.bool)
    no_ifgs, rows, cols = nan_ifg.shape

    for r in range(rows):
        for c in range(cols):
            nan_v = nan_ifg[ifg_list_mst_id, r, c]
            # if there is nan value in the independent ifglist, redo mst search
            if np.count_nonzero(nan_v) > 0:
                nan_v = nan_ifg[:, r, c]
                if no_ifgs - np.count_nonzero(nan_v) >= p_thresh_hold:
                    # get all valid ifgs from ifglist,
                    # and then select the ones that are not nan on this pixel
                    ifg_list_valid = get_sub_structure(ifg_list, nan_v)
                    # calculate mst again
                    ifglist_mst_valid_id = matlab_mst_kruskal(ifg_list_valid)
                    mst_mat[ifglist_mst_valid_id, r, c] = 1
            else:
                mst_mat[ifg_list_mst_id, r, c] = 1
    return mst_mat


def get_sub_structure(ifg_list, nan_v):
    """
    This is the getsucstruct.m in pi-rate/matlab.
    :param ifg_list: original ifg_list class instance.
    :param nan_v: all ifg values at this location.
    :return:
    """
    # TODO: Matlab does not pass through get_nml a second time.
    # Check which is correct. I think running via get_nml here is safer
    data_files_valid = [ifg_list.nml[a] for a in np.nonzero(~nan_v)[0]]
    ifg_list_valid, _ = get_nml(datafiles=data_files_valid,
                                nan_conversion=True)
    return ifg_list_valid

if __name__ == "__main__":
    ifg_list, epoch_list = get_nml(nan_conversion=True)
    ifg_list_mst_id = matlab_mst_kruskal(ifg_list)
    matlab_mst(ifg_list, ifg_list_mst_id)
