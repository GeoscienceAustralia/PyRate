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
This Python module implements the minimum spanning tree
functionality for selecting interferometric observations.
"""
# pylint: disable=invalid-name
from itertools import product
from numpy import array, nan, isnan, float32, empty, sum as nsum
import numpy as np
from networkx.classes.reportviews import EdgeView
import networkx as nx
from joblib import Parallel, delayed

from pyrate.core.algorithm import ifg_date_lookup
from pyrate.core.algorithm import ifg_date_index_lookup
from pyrate.core import config as cf
from pyrate.core.shared import IfgPart, create_tiles
from pyrate.core.shared import joblib_log_level
from pyrate.core.logger import pyratelogger as log

np.seterr(invalid='ignore')  # stops RuntimeWarning in nan conversion

# TODO: may need to implement memory saving row-by-row access
# TODO: document weighting by either Nan fraction OR variance




def mst_from_ifgs(ifgs):
    """
    Returns a Minimum Spanning Tree (MST) network from the given interferograms
    using Kruskal's algorithm. The MST is calculated using a weighting based
    on the number of NaN cells in the phase band.

    :param list ifgs: List of interferogram objects (Ifg class)

    :return: edges: The number of network connections
    :rtype: int
    :return: is_tree: A boolean that is True if network is a tree
    :rtype: bool
    :return: ntrees: The number of disconnected trees
    :rtype: int
    :return: mst_ifgs: Minimum Spanning Tree network of interferograms
    :rtype: list
    """

    edges_with_weights_for_networkx = [(i.master, i.slave, i.nan_fraction)
                                       for i in ifgs]
    g_nx = _build_graph_networkx(edges_with_weights_for_networkx)
    mst = nx.minimum_spanning_tree(g_nx)
    # mst_edges, is tree?, number of trees
    edges = mst.edges()
    ifg_sub = [ifg_date_index_lookup(ifgs, d) for d in edges]
    mst_ifgs = [i for k, i in enumerate(ifgs) if k in ifg_sub]
    return mst.edges(), nx.is_tree(mst), \
        nx.number_connected_components(mst), mst_ifgs


def mst_parallel(ifgs, params):
    """
    Wrapper function for calculating MST matrix in non-MPI runs.

    :param list ifgs: List of interferogram objects (Ifg class)
    :param dict params: Dictionary of parameters

    :return: result: 3-dimensional Minimum Spanning Tree matrix
    :rtype: ndarray
    """

    log.info('Calculating MST in tiles')
    ncpus = params[cf.PROCESSES]
    no_ifgs = len(ifgs)
    no_y, no_x = ifgs[0].phase_data.shape
    tiles = create_tiles(ifgs[0].shape)
    no_tiles = len(tiles)
    # need to break up the ifg class as multiprocessing does not allow pickling
    # don't read in all the phase data at once
    # use data paths and locally read in each core/process,
    # gdal read is very fast
    ifg_paths = [i.data_path for i in ifgs]
    result = empty(shape=(no_ifgs, no_y, no_x), dtype=np.bool)

    if params[cf.PARALLEL]:
        log.info('Calculating MST using {} tiles in parallel using {} ' \
                 'processes'.format(no_tiles, ncpus))
        t_msts = Parallel(n_jobs=params[cf.PROCESSES], 
                          verbose=joblib_log_level(cf.LOG_LEVEL))(
            delayed(mst_multiprocessing)(t, ifg_paths, params=params)
            for t in tiles)
        for k, tile in enumerate(tiles):
            result[:, tile.top_left_y:tile.bottom_right_y,
                   tile.top_left_x: tile.bottom_right_x] = t_msts[k]
    else:
        log.info('Calculating MST using {} tiles in serial'.format(no_tiles))
        for k, tile in enumerate(tiles):
            result[:, tile.top_left_y:tile.bottom_right_y,
                   tile.top_left_x: tile.bottom_right_x] = \
                mst_multiprocessing(tile, ifg_paths, params=params)

    return result


def mst_multiprocessing(tile, ifgs_or_paths, preread_ifgs=None, params=None):
    """
    Wrapper function for calculating MST matrix for a tile

    :param ndarray tile: Tile class instance
    :param list ifgs_or_paths: All interferograms paths of the problem.
        List of strings
    :param dict preread_ifgs: Dictionary of interferogram metadata

    :return: mst_tile: MST matrix tile. An array of booleans representing
        valid ifg connections
    :rtype: ndarray
    """
    #The memory requirement during MPI MST computation is determined by the
    #number of interferograms times size of IfgPart. Note that we need all
    #interferogram header information (like masters/slave dates) for MST
    #computation. To manage memory we need smaller tiles (IfgPart) as number
    #of interferograms increases

    ifg_parts = [IfgPart(p, tile, preread_ifgs, params) for p in ifgs_or_paths]
    return mst_boolean_array(ifg_parts)


def _build_graph_networkx(edges_with_weights):
    """
    Convenience graph builder function: returns a new graph object.
    """
    g = nx.Graph()
    g.add_weighted_edges_from(edges_with_weights)
    return g


def mst_boolean_array(ifgs):
    """
    Returns a 3D array of booleans constituting valid interferogram connections
    in the Minimum Spanning Tree matrix.

    :param list ifgs: Sequence of interferogram objects

    :return: result: Array of booleans representing valid ifg connections
    :rtype: ndarray
    """
    #The MSTs are stripped of connecting edge info, leaving just the ifgs.
    nifgs = len(ifgs)
    ny, nx = ifgs[0].phase_data.shape
    result = empty(shape=(nifgs, ny, nx), dtype=np.bool)
    
    for y, x, mst in mst_matrix_networkx(ifgs):
        # mst is a list of datetime.date tuples
        if isinstance(mst, EdgeView):
            ifg_sub = [ifg_date_index_lookup(ifgs, d) for d in mst]
            ifg_sub_bool = [True if i in ifg_sub else False
                            for i in range(nifgs)]  # boolean conversion
            result[:, y, x] = np.array(ifg_sub_bool)
        else:
            result[:, y, x] = np.zeros(nifgs, dtype=bool)
    return result


def _mst_matrix_ifgs_only(ifgs):
    """
    Alternative method for producing 3D MST array
    """
    #Currently not used
    #The MSTs are stripped of connecting edge info, leaving just the ifgs.
    result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix_networkx(ifgs):
        if isinstance(mst, EdgeView):
            ifg_sub = [ifg_date_lookup(ifgs, d) for d in mst]
            result[(y, x)] = tuple(ifg_sub)
        else:
            result[(y, x)] = mst  # usually NaN
    return result


def _mst_matrix_as_array(ifgs):
    """
    Alternative method for producing 3D MST array
    """
    #Currently not used
    #Each pixel contains an MST (with connecting edges etc).
    mst_result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix_networkx(ifgs):
        mst_result[y, x] = mst
    return mst_result


# TODO: custom edge weighting could be included with an additional
# 'weights' arg if some other weighting criterion is required later
def mst_matrix_networkx(ifgs):
    """
    Generates MST network for a single pixel for the given ifgs using
    NetworkX-package algorithms.

    :param list ifgs: Sequence of interferogram objects

    :return: y: pixel y coordinate
    :rtype: int
    :return: x: pixel x coordinate
    :rtype: int
    :return: mst: list of tuples for edges in the minimum spanning tree
    :rtype: list
    """
    # make default MST to optimise result when no Ifg cells in a stack are nans
    edges_with_weights = [(i.master, i.slave, i.nan_fraction) for i in ifgs]
    edges, g_nx = _minimum_spanning_edges_from_mst(edges_with_weights)
    # TODO: memory efficiencies can be achieved here with tiling

    list_of_phase_data = [i.phase_data for i in ifgs]
    log.debug("list_of_phase_data length: " + str(len(list_of_phase_data)))
    for row in list_of_phase_data:
        log.debug("row length in list_of_phase_data: " + str(len(row)))
        log.debug("row in list_of_phase_data: " + str(row))
    data_stack = array(list_of_phase_data, dtype=float32)

    # create MSTs for each pixel in the ifg data stack
    nifgs = len(ifgs)

    for y, x in product(range(ifgs[0].nrows), range(ifgs[0].ncols)):
        values = data_stack[:, y, x]  # vertical stack of ifg values for a pixel
        nan_count = nsum(isnan(values))

        # optimisations: use pre-created results for all nans/no nans
        if nan_count == 0:
            yield y, x, edges
            continue
        elif nan_count == nifgs:
            yield y, x, nan
            continue

        # dynamically modify graph to reuse a single graph: this should avoid
        # repeatedly creating new graph objs & reduce RAM use
        ebunch_add = []
        ebunch_delete = []
        for value, edge in zip(values, edges_with_weights):
            if not isnan(value):
                if not g_nx.has_edge(edge[0], edge[1]):
                    ebunch_add.append(edge)
            else:
                if g_nx.has_edge(edge[0], edge[1]):
                    ebunch_delete.append(edge)
        if ebunch_add:
            g_nx.add_weighted_edges_from(ebunch_add)
        if ebunch_delete:
            g_nx.remove_edges_from(ebunch_delete)
        yield y, x, nx.minimum_spanning_tree(g_nx).edges()


def _minimum_spanning_edges_from_mst(edges):
    """
    Convenience function to determine MST edges
    """
    g_nx = _build_graph_networkx(edges)
    T = nx.minimum_spanning_tree(g_nx)  # step ifglist_mst in make_mstmat.m
    edges = T.edges()
    return edges, g_nx
