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
This Python module implements the minimum spanning tree
functionality for selecting interferometric observations.
"""
# pylint: disable=invalid-name
from __future__ import print_function
import logging
from itertools import product
from numpy import array, nan, isnan, float32, empty, sum as nsum
import numpy as np
import networkx as nx
from joblib import Parallel, delayed

from pyrate.algorithm import ifg_date_lookup
from pyrate.algorithm import ifg_date_index_lookup
from pyrate import config as cf
from pyrate.shared import IfgPart, create_tiles
np.seterr(invalid='ignore')  # stops RuntimeWarning in nan conversion

# TODO: may need to implement memory saving row-by-row access
# TODO: document weighting by either Nan fraction OR variance

log = logging.getLogger(__name__)


def mst_from_ifgs(ifgs):
    """
    Returns default MST dict for the given Ifgs. The MST is calculated using a
    weighting based on the number of incoherent cells in the phase band.

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
    wrapper function for calculating ifgs in non mpi runs
    :param ifgs:
    :param params:
    :return:
    """
    print('Calculating mst using tiles')
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
        print('Calculating mst using {} tiles in parallel using {} ' \
              'processes'.format(no_tiles, ncpus))
        t_msts = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(mst_multiprocessing)(t, ifg_paths)
            for t in tiles)
        for k, tile in enumerate(tiles):
            result[:, tile.top_left_y:tile.bottom_right_y,
                   tile.top_left_x: tile.bottom_right_x] = t_msts[k]
    else:
        print('Calculating mst using {} tiles in serial'.format(no_tiles))
        for k, tile in enumerate(tiles):
            result[:, tile.top_left_y:tile.bottom_right_y,
                   tile.top_left_x: tile.bottom_right_x] = \
                mst_multiprocessing(tile, ifg_paths)

    return result


def mst_multiprocessing(tile, ifgs_or_paths, preread_ifgs=None):
    """
    The memory requirement during mpi mst computation is determined by the
    number of ifgs times size of IfgPart. Note that we need all ifg header
    information (like masters/slave dates) for mst computation.
    To manage memory we need smaller tiles (IfgPart) as number of ifgs go up.
    :param tile: Tile class instance
    :param ifgs_or_paths: all ifg paths of the problem. List of strings.
    """
    ifg_parts = [IfgPart(p, tile, preread_ifgs) for p in ifgs_or_paths]
    t_mst = mst_boolean_array(ifg_parts)
    return t_mst


def _build_graph_networkx(edges_with_weights):
    """
    Convenience graph builder function: returns a new graph obj.
    """
    g = nx.Graph()
    g.add_weighted_edges_from(edges_with_weights)
    return g


def mst_boolean_array(ifgs):
    """
    Filter: returns array of independent ifgs from the pixel by pixel MST,
    like that used by the MATLAB version Pirate.

    The MSTs are stripped of connecting edge info, leaving just the ifgs.

    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    """
    no_ifgs = len(ifgs)
    no_y, no_x = ifgs[0].phase_data.shape
    result = empty(shape=(no_ifgs, no_y, no_x), dtype=np.bool)

    for y, x, mst in mst_matrix_networkx(ifgs):
        # mst is a list of datetime.date tuples
        if isinstance(mst, list):
            ifg_sub = [ifg_date_index_lookup(ifgs, d) for d in mst]
            ifg_sub_bool = [True if i in ifg_sub else False
                            for i in range(no_ifgs)]  # boolean conversion
            result[:, y, x] = np.array(ifg_sub_bool)
        else:
            result[:, y, x] = np.zeros(no_ifgs, dtype=bool)
    return result


def mst_matrix_ifgs_only(ifgs):
    """
    Filter: returns array of independent ifgs from the pixel by pixel MST.

    The MSTs are stripped of connecting edge info, leaving just the ifgs.

    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    """

    result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix_networkx(ifgs):
        if isinstance(mst, list):
            ifg_sub = [ifg_date_lookup(ifgs, d) for d in mst]
            result[(y, x)] = tuple(ifg_sub)
        else:
            result[(y, x)] = mst  # usually NaN
    return result


def mst_matrix_as_array(ifgs):
    """
    Filter: returns array of pixel by pixel MSTs.

    Each pixel contains an MST (with connecting edges etc).

    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    """

    mst_result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix_networkx(ifgs):
        mst_result[y, x] = mst
    return mst_result


# TODO: custom weighting could included with an additional 'weights' arg if some
# other weighting criterion is required later
def mst_matrix_networkx(ifgs):
    """
    Generates/emits MST trees on a pixel-by-pixel basis for the given ifgs.
    :param ifgs: sequence of Ifg objs
    """

    # make default MST to optimise result when no Ifg cells in a stack are nans
    edges_with_weights = [(i.master, i.slave, i.nan_fraction) for i in ifgs]
    edges, g_nx = minimum_spanning_edges_from_mst(edges_with_weights)

    # TODO: memory efficiencies can be achieved here with tiling
    data_stack = array([i.phase_data for i in ifgs], dtype=float32)

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
        yield y, x, minimum_spanning_tree(g_nx).edges()


def minimum_spanning_edges_from_mst(edges):
    """
    Parameters
    ----------
    edges: list
        List of tuples (master, slave, nan_frac) corresponding to ifgs
    Returns
    -------
    edges: list
        list of mst edges
    g_nx: nx.Graph() instance
    """
    g_nx = _build_graph_networkx(edges)
    T = minimum_spanning_tree(g_nx)  # step ifglist_mst in make_mstmat.m
    edges = T.edges()
    return edges, g_nx


def minimum_spanning_tree(G, weight='weight'):
    """Return a minimum spanning tree or forest of an undirected
    weighted graph.

    A minimum spanning tree is a subgraph of the graph (a tree) with
    the minimum sum of edge weights.

    If the graph is not connected a spanning forest is constructed.  A
    spanning forest is a union of the spanning trees for each
    connected component of the graph.

    Parameters
    ----------
    G : NetworkX Graph

    weight : string
       Edge data key to use for weight (default 'weight').

    Returns
    -------
    G : NetworkX Graph
       A minimum spanning tree or forest.

    Examples
    --------
    >>> G=nx.cycle_graph(4)
    >>> G.add_edge(0,3,weight=2) # assign weight 2 to edge 0-3
    >>> T=nx.minimum_spanning_tree(G)
    >>> print(sorted(T.edges(data=True)))
    [(0, 1, {}), (1, 2, {}), (2, 3, {})]

    Notes
    -----
    Uses Kruskal's algorithm.

    If the graph edges do not have a weight attribute a default weight of 1
    will be used.
    """
    T = nx.Graph(minimum_spanning_edges(G, weight=weight, data=True))
    # Add isolated nodes
    if len(T) != len(G):
        T.add_nodes_from([n for n, d in G.degree().items() if d == 0])
    # Add node and graph attributes as shallow copy
    for n in T:
        T.node[n] = G.node[n].copy()
    T.graph = G.graph.copy()
    return T


def minimum_spanning_edges(G, weight='weight', data=True):
    """Generate edges in a minimum spanning forest of an undirected
    weighted graph.

    A minimum spanning tree is a subgraph of the graph (a tree)
    with the minimum sum of edge weights.  A spanning forest is a
    union of the spanning trees for each connected component of the graph.

    Parameters
    ----------
    G : NetworkX Graph

    weight : string
       Edge data key to use for weight (default 'weight').

    data : bool, optional
       If True yield the edge data along with the edge.

    Returns
    -------
    edges : iterator
       A generator that produces edges in the minimum spanning tree.
       The edges are three-tuples (u,v,w) where w is the weight.

    Examples
    --------
    >>> G=nx.cycle_graph(4)
    >>> G.add_edge(0,3,weight=2) # assign weight 2 to edge 0-3
    >>> mst=nx.minimum_spanning_edges(G,data=False) # a generator of MST edges
    >>> edgelist=list(mst) # make a list of the edges
    >>> print(sorted(edgelist))
    [(0, 1), (1, 2), (2, 3)]

    Notes
    -----
    Uses Kruskal's algorithm.

    If the graph edges do not have a weight attribute a default weight of 1
    will be used.

    Modified code from David Eppstein, April 2006
    http://www.ics.uci.edu/~eppstein/PADS/
    """
    # Modified code from David Eppstein, April 2006
    # http://www.ics.uci.edu/~eppstein/PADS/
    # Kruskal's algorithm: sort edges by weight, and add them one at a time.
    # We use Kruskal's algorithm, first because it is very simple to
    # implement once UnionFind exists, and second, because the only slow
    # part (the sort) is sped up by being built in to Python.
    from networkx.utils import UnionFind
    if G.is_directed():
        raise nx.NetworkXError(
            "Mimimum spanning tree not defined for directed graphs.")

    subtrees = UnionFind()
    edges = sorted(G.edges(data=True), key=lambda t: t[2].get(weight, 1.0))

    for u, v, d in edges:
        if subtrees[u] != subtrees[v]:
            if data:
                yield (u, v, d)
            else:
                yield (u, v)
            subtrees.union(u, v)
