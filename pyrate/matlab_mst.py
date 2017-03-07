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
This Python module implements the minimum spanning tree matrix algorithm
of the Matlab code in 'make_mstmat.m' of the Pirate package.
"""
# pylint: disable=missing-docstring
from __future__ import print_function
import itertools
from abc import ABCMeta, abstractmethod

import numpy as np

from pyrate.algorithm import get_epochs
from pyrate.shared import Ifg

DTYPE = [('id', int), ('master', int), ('slave', int), ('nan_frac', float)]


class IfGMeta(object):
    """
    Metaclass for matlab mst calculation
    """
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

    def update_nan_frac(self, nodata):
        for i in self.ifgs:
            i.nodata_value = nodata
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
    """
    copy of matlab ifglist in getnml.m.
    Please note that we don't need BaseT unless we are using variance in ifg
    data as cost.
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, datafiles=None):
        self.datafiles = datafiles
        self.nml = self.get_nml_list()
        self.ifgs = self.get_ifgs_list()
        self.id = range(len(self.nml))
        self.base_t = None
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


def data_setup(datafiles):
    '''Returns Ifg objs for the files in the sydney test dir
    input phase data is in radians;
    these ifgs are in radians - not converted to mm'''
    datafiles.sort()
    ifgs = [Ifg(i) for i in datafiles]

    for i in ifgs:
        i.open()
    return ifgs


def get_nml(ifg_list_instance, nodata_value,
            nan_conversion=False):
    """
    A reproduction of getnml.m, the matlab function in pi-rate.
    Note: the matlab version tested does not have nan's.
    replaces the ifg_list_instance in place
    """
    _epoch_list, n = get_epochs(ifg_list_instance.ifgs)
    ifg_list_instance.reshape_n(n)
    if nan_conversion:
        ifg_list_instance.update_nan_frac(nodata_value)
        # turn on for nan conversion
        ifg_list_instance.convert_nans(nan_conversion=nan_conversion)
    ifg_list_instance.make_data_stack()
    return ifg_list_instance, _epoch_list


def matlab_mst_kruskal(edges, ntrees=False):
    """
    This is an implementation of the pi-rate mst_kruskal.m
    :param edges: list of edges, list of tuples (id, master, slave, nan_frac)
    :return:
    """

    num_ifgs = len(edges)
    master_l = [e[1] for e in edges]
    slave_l = [e[2] for e in edges]

    num_images = max(max(master_l), max(slave_l))

    # sort edges based on nan_frac
    ifg_sorted = sorted(edges, key=lambda t: t[3])

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
    if ntrees:
        return calculate_connect_and_ntrees(connect, mst_list)
    else:
        return [i[0] for i in mst_list]


def calculate_connect_and_ntrees(connect, mst_list):
    """
    Count isolated trees
    """
    zero_count = np.where(np.sum(connect, axis=1) == 0)[0]
    if zero_count.shape[0]:
        raise ValueError('Input connect matrix not compatible')
    cnt = np.where(np.sum(connect, axis=1) == 1)
    connect = np.delete(connect, cnt, axis=0)
    ntrees = connect.shape[0]
    if not ntrees:
        raise ValueError('No tree found. Check connect matrix')
    # return mst, connect, ntrees
    return [i[0] for i in mst_list], connect, ntrees


def matlab_mst(ifg_object, p_threshold=1):
    """
    This is an implementation of matlab/pirate make_mstmat.m.
    """
    edges = get_sub_structure(ifg_object,
                              np.zeros(len(ifg_object.id), dtype=bool))

    ifg_list_mst_id = matlab_mst_kruskal(edges)
    data_stack = ifg_object.data_stack
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
                    edges = get_sub_structure(ifg_object, nan_v)
                    # calculate mst again
                    ifglist_mst_valid_id = matlab_mst_kruskal(edges)
                    mst_mat[ifglist_mst_valid_id, r, c] = 1
                else:
                    pass
            else:
                mst_mat[ifg_list_mst_id, r, c] = 1
    return mst_mat


def matlab_mst_gen(ifg_instance, p_threshold=1):

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
    edges = get_sub_structure(ifg_instance,
                              np.zeros(len(ifg_instance.id), dtype=bool))

    ifg_list_mst_id = matlab_mst_kruskal(edges)
    data_stack = ifg_instance.data_stack
    # np.array([i.phase_data for i in ifg_list.ifgs],
    #                   dtype=np.float32)
    # nan_ifg = np.isnan(data_stack)  # doing nan conversion later saves memory
    num_ifgs, rows, cols = data_stack.shape

    for r, c in itertools.product(range(rows), range(cols)):
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
                new_edges = get_sub_structure(
                    ifg_instance, nan_v)
                # calculate mst again
                ifglist_mst_valid_id = matlab_mst_kruskal(new_edges)
                mst_yield[ifglist_mst_valid_id] = True
                yield r, c, mst_yield
            else:
                yield r, c, mst_yield


def matlab_mst_bool(ifg_list_instance, p_threshold=1):
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

    for y, x, mst in matlab_mst_gen(ifg_list_instance,
                                    p_threshold):
        result[:, y, x] = mst
    return result


def get_sub_structure(ifg_list, nan_v):
    """
    This is the getsucstruct.m in pi-rate/matlab.
    :param ifg_list: original ifg_list class instance.
    :param nan_v: all ifg values at this location.
    Returns list of tuples (id, master, slave, nan_frac)
    corresponding to the chosen ifgs.
    """
    indices_chosen = np.nonzero(~nan_v)[0]

    return [(ifg_list.id[i],
             ifg_list.master_num[i],
             ifg_list.slave_num[i],
             ifg_list.nan_frac[i])
            for i in indices_chosen]
