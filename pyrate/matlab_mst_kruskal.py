__author__ = 'Sudipta Basak'
__date_created__ = '7/12/15'

"""
This is the python implementation of the make_mstmat.m very closely
resembling the matlab.
"""

import sys
import numpy as np
import itertools
from abc import ABCMeta, abstractmethod

from pyrate.tests.common import sydney_data_setup_ifg_file_list
from pyrate.tests.common import sydney_data_setup
from pyrate.algorithm import get_epochs_and_n
from pyrate.shared import Ifg

DTYPE = [('id', int), ('master', int), ('slave', int), ('nan_frac', float)]

class IfGMeta(object):
    __metaclass__ = ABCMeta
    ifgs = None
    id = None
    data_stack = None
    master_num = None
    slave_num = None
    n = None
    nan_frac = None

    @abstractmethod
    def get_nml_list(self):
        """
        Override this method to get the data files.
        """
        return

    @abstractmethod
    def get_ifgs_list(self):
        """
        Override this method to get the difgs from the datafiles.
        """
        return

    def reshape_n(self, n):
        n_reshaped = np.reshape(n, newshape=(2, len(self.id)))
        self.master_num = n_reshaped[0, :]
        self.slave_num = n_reshaped[1, :]
        self.n = n

    def update_nan_frac(self):
        self.nan_frac = [i.nan_fraction for i in self.ifgs]

    def convert_nans(self, nan_conversion=False):
        if nan_conversion:
            for i in self.ifgs:
                if not i.nan_converted:
                    i.convert_to_nans()

    def make_data_stack(self):
        self.data_stack = np.array([i.phase_data for i in self.ifgs],
                                   dtype=np.float32)


class IfgListPyRate(IfGMeta):

    def __init__(self, datafiles):
        self.datafiles = datafiles
        self.nml = self.get_nml_list()
        self.ifgs = self.get_ifgs_list()
        self.id = range(len(self.nml))
        self.base_T = None
        self.max_var = np.zeros_like(self.id)
        self.alpha = np.zeros_like(self.id)
        self.nan_frac = np.zeros_like(self.id)
        self.master_num = None
        self.slave_num = None
        self.n = None

    def get_ifgs_list(self):
        return data_setup(self.datafiles)

    def get_nml_list(self):
        return self.datafiles


class IfgListMatlabTest(IfGMeta):
    """
    copy of matlab ifglist in getnml.m.
    1. Please note that we don't need BaseT unless we are using variance in ifg
    data as cost.
    2.
    """

    def __init__(self, datafiles=None):
        self.datafiles = datafiles
        self.nml = self.get_nml_list()
        self.ifgs = self.get_ifgs_list()
        self.id = range(len(self.nml))
        self.base_T = None
        self.max_var = np.zeros_like(self.id)
        self.alpha = np.zeros_like(self.id)
        self.nan_frac = np.zeros_like(self.id)
        self.master_num = None
        self.slave_num = None
        self.n = None

    def get_ifgs_list(self):
        return sydney_data_setup(self.datafiles)

    def get_nml_list(self):
        return sydney_data_setup_ifg_file_list(self.datafiles)


def data_setup(datafiles):
    '''Returns Ifg objs for the files in the sydney test dir
    input phase data is in radians;
    these ifgs are in radians - not converted to mm'''
    datafiles.sort()
    ifgs = [Ifg(i) for i in datafiles]

    for i in ifgs:
        i.open()

    return ifgs


def get_nml(ifg_list_instance, nan_conversion=False, prefix_len=4):
    """
    A reproduction of getnml.m, the matlab function in pi-rate.
    Note: the matlab version tested does not have nan's.
    """
    epoch_list_, n = get_epochs_and_n(ifg_list_instance.ifgs)
    ifg_list_instance.reshape_n(n)
    if nan_conversion:
        ifg_list_instance.update_nan_frac()  # turn on for nan conversion
        ifg_list_instance.convert_nans(nan_conversion=nan_conversion)
    ifg_list_instance.make_data_stack()
    return ifg_list_instance, epoch_list_


def sort_list(id_l, master_l, slave_l, nan_frac_l):
    sort_list = map(lambda i, m, s, n: (i, m, s, n),
                    id_l, master_l, slave_l, nan_frac_l)

    sort_list = np.array(sort_list, dtype=DTYPE)
    return np.sort(sort_list, order=['nan_frac'])


def matlab_mst_kruskal(id_l, master_l, slave_l, nan_frac_l, connect_flag=False):
    """
    This is an implementation of the pi-rate mst_kruskal.m
    :param id_l: list of ifg file ids
    :param master_l: list of ifg master dates
    :param slave_l: list of ifg slave dates
    :param nan_frac_l: list of ifg nan fractions
    :return:
    """

    num_ifgs = len(master_l)
    num_images = max(max(master_l), max(slave_l))
    # print 'ifg_list', ifg_list_
    ifg_sorted = sort_list(id_l, master_l, slave_l, nan_frac_l)
    # print 'ifg_sorted', ifg_sorted

    # add one to ensure index number + 1
    connect = np.eye(num_images + 1, dtype=np.bool)

    mst_list = []

    for i in range(num_ifgs):
        master = ifg_sorted[i][1]
        slave = ifg_sorted[i][2]
        loc_master = np.where(connect[:, master] == 1)[0][0]
        loc_slave = np.where(connect[:, slave] == 1)[0][0]

        if loc_master != loc_slave:
            mst_list.append(ifg_sorted[i])
            connect[loc_master, :] = connect[loc_master, :] + \
                                     connect[loc_slave, :]
            connect = np.delete(connect, loc_slave, axis=0)

    mst_list = np.array(mst_list, dtype=DTYPE)
    mst_list = np.sort(mst_list, order=['id'])

    # count isolated trees
    # TODO: test connect_flag block below
    if connect_flag:
        cnt = np.where(np.sum(connect, axis=1) == 1)
        connect = np.delete(connect, cnt, axis=0)
        # return mst, connect, ntrees
        return [i[0] for i in mst_list], connect, connect.shape[0]
    else:
        return [i[0] for i in mst_list]


def matlab_mst(ifg_list_, p_threshold=1):
    """
    This is an implementation of matlab/pirate make_mstmat.m.
    """
    ifg_list_mst_id = matlab_mst_kruskal(ifg_list_.id,
        ifg_list_.master_num, ifg_list_.slave_num, ifg_list_.nan_frac)
    data_stack = ifg_list_.data_stack
    nan_ifg = np.isnan(data_stack)
    mst_mat = np.zeros_like(nan_ifg, dtype=np.bool)
    num_ifgs, rows, cols = nan_ifg.shape

    for r in range(rows):
        for c in range(cols):
            nan_v = nan_ifg[ifg_list_mst_id, r, c]
            # if there is nan value in the independent ifglist, redo mst search
            if np.count_nonzero(nan_v) > 0:
                nan_v = nan_ifg[:, r, c]
                if num_ifgs - np.count_nonzero(nan_v) >= p_threshold:
                    # get all valid ifgs from ifglist,
                    # and then select the ones that are not nan on this pixel
                    id, master, slave, nan_frac = get_sub_structure(
                        ifg_list_, nan_v)
                    # calculate mst again
                    ifglist_mst_valid_id = matlab_mst_kruskal(id, master,
                                                              slave, nan_frac)
                    mst_mat[ifglist_mst_valid_id, r, c] = 1
                else:
                    pass
            else:
                mst_mat[ifg_list_mst_id, r, c] = 1
    return mst_mat


def matlab_mst_generator_boolean_array(ifg_instance, p_threshold=1):

    """
    :param ifg_instance: IfgListPyRate instance
    :param p_threshold: minimum number of non-nan values at any pixel for selection

    This is an implementation of matlab/pirate make_mstmat.m.
    This we will be able to call from mst.py and the rest of the
    python framework setup so far.

    Please note that the generator version is more memory efficient.
    If memory was not a concern we could have found the entire mst matrix in the
    previous function and this would have been unnecessary.
    """
    ifg_list_mst_id = matlab_mst_kruskal(ifg_instance.id,
        ifg_instance.master_num, ifg_instance.slave_num, ifg_instance.nan_frac)
    data_stack = ifg_instance.data_stack
    # np.array([i.phase_data for i in ifg_list.ifgs],
    #                   dtype=np.float32)
    # nan_ifg = np.isnan(data_stack)  # doing nan conversion later saves memory
    num_ifgs, rows, cols = data_stack.shape

    for r, c in itertools.product(xrange(rows), xrange(cols)):
        # nan_count of the vertical stack of ifg values for a pixel
        mst_yield = np.zeros(num_ifgs, dtype=np.bool)
        nan_count = np.sum(np.isnan(data_stack[ifg_list_mst_id, r, c]))
        # if there is nan value in the independent ifglist, redo mst search
        if nan_count == 0:
            mst_yield[ifg_list_mst_id] = True
            yield r, c, mst_yield
        else:
            nan_v = np.isnan(data_stack[:, r, c])
            nan_count = np.sum(nan_v)            
            if (num_ifgs - nan_count) >= p_threshold:
                # get all valid ifgs from ifglist,
                # and then select the ones that are not nan on this pixel
                id, master, slave, nan_frac = get_sub_structure(
                    ifg_instance, nan_v)
                # calculate mst again
                ifglist_mst_valid_id = matlab_mst_kruskal(
                    id, master, slave, nan_frac)
                mst_yield[ifglist_mst_valid_id] = True
                yield r, c, mst_yield
            else:
                yield r, c, mst_yield

# TODO: performance test matlab_mst_boolean_array vs matlab_mst for large ifgs
def matlab_mst_boolean_array(ifg_list_instance, p_threshold=1):
    """
    :param ifg_instance: IfgListPyRate instance
    :param p_threshold: minimum number of non-nan values at any pixel for selection

    This should have the same output as matlab_mst. Should be tested.
    Please note that the generator version is more memory efficient.
    If memory was not a concern we could have found the entire mst matrix in the
    previous function and this would have been unnecessary.
    :return:
    """
    num_ifgs = len(ifg_list_instance.ifgs)
    no_y, no_x = ifg_list_instance.ifgs[0].phase_data.shape
    result = np.empty(shape=(num_ifgs, no_y, no_x), dtype=np.bool)

    for y, x, mst in matlab_mst_generator_boolean_array(ifg_list_instance,
                                                        p_threshold):
        result[:, y, x] = mst
    return result


def matlab_mst_kruskal_from_ifgs(ifgs):
    dest_paths = [i.data_path for i in ifgs]
    ifg_instance = IfgListPyRate(datafiles=dest_paths)
    ifg_instance_updated, epoch_list = \
        get_nml(ifg_instance, nan_conversion=True)

    ifg_list_mst_id = matlab_mst_kruskal(
        ifg_instance_updated.id, ifg_instance_updated.master_num,
        ifg_instance_updated.slave_num, ifg_instance_updated.nan_frac)

    return [ifgs[i] for i in ifg_list_mst_id]


def str_to_class(str):
    """
    :param str: a string
    :return: class from string
    """
    return getattr(sys.modules[__name__], str)


def get_all_class_attributes(this_class):
    print dir(this_class)
    return [attr for attr in dir(this_class)
            if not callable(attr)
            and not attr.startswith("_")]


def get_all_attriblues_of_class(class_instance):
    return vars(class_instance).keys()


def get_sub_structure(ifg_list, nan_v):
    """
    This is the getsucstruct.m in pi-rate/matlab.
    :param ifg_list: original ifg_list class instance.
    :param nan_v: all ifg values at this location.
    :return:
    """
    indices_chosen = np.nonzero(~nan_v)[0]

    # TODO: remove the list comprehensions
    id = [ifg_list.id[i] for i in indices_chosen]
    master_num = [ifg_list.master_num[i] for i in indices_chosen]
    slave_num = [ifg_list.slave_num[i] for i in indices_chosen]
    nan_frac = [ifg_list.nan_frac[i] for i in indices_chosen]

    return id, master_num, slave_num, nan_frac

if __name__ == "__main__":
    ifg_instance_main = IfgListMatlabTest()
    ifg_list, epoch_list = get_nml(ifg_instance_main, nan_conversion=True)
    mst_mat1 = matlab_mst(ifg_list)
    mst_mat2 = matlab_mst_boolean_array(ifg_list)
    print np.array_equal(mst_mat1, mst_mat2)  # assert equality of both methods
