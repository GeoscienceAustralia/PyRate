'''
Minimum Spanning Tree functionality for PyRate.

Contains functions to calculate MST using interferograms.

.. codeauthor:: Ben Davies
'''

from itertools import product
from numpy import array, nan, isnan, float32, empty

from pyrate.algorithm import get_all_epochs, master_slave_ids, ifg_date_lookup, get_epochs
from pyrate.algorithm import ifg_date_index_lookup

from pygraph.classes.graph import graph
from pygraph.algorithms.minmax import minimal_spanning_tree
import networkx as nx

# TODO: may need to implement memory saving row-by-row access
# TODO: document weighting by either Nan fraction OR variance


def _remove_root_node(mst):
    """
    Discard pygraph's root nodes from MST dict to conserve memory.
    Note multiple root nodes exist in MST dict.
    """
    for k in mst.keys():
        if mst[k] is None:
            del mst[k]
    return mst


def _remove_root_node_new(mst):
    """
    Discard pygraph's root nodes from MST dict to conserve memory.
    Note multiple root nodes exist in MST dict.
    """
    return {k: v for k, v in mst.items() if v}


# TODO: this does not appear to be used anywhere
def default_mst(ifgs, noroot=True):
    '''
    Returns default MST dict for the given Ifgs. The MST is calculated using a
    weighting based on the number of incoherent cells in the phase band.

    :param noroot: True removes the PyGraph default root node from the result.
    '''

    edges = [(i.master, i.slave) for i in ifgs]
    epochs = master_slave_ids(get_all_epochs(ifgs)).keys()
    weights = [i.nan_fraction for i in ifgs]  # NB: other attrs for weights?

    g = _build_graph(epochs, edges, weights)
    mst = minimal_spanning_tree(g)
    return _remove_root_node(mst) if noroot else mst


def _build_graph(nodes, edges, weights, noroot=True):
    """
    Convenience graph builder function: returns a new graph obj.
    """

    g = graph()
    g.add_nodes(nodes)
    for edge, weight in zip(edges, weights):
        g.add_edge(edge, wt=weight)
    return g


def _build_graph_networkx(nodes, edges_with_weights, noroot=True):
    """
    Convenience graph builder function: returns a new graph obj.
    """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_weighted_edges_from(edges_with_weights)
    return g


def _build_graph_new(nodes, edges, weights, noroot=True):
    """
    Convenience graph builder function: returns a new graph obj.
    Equivalent performance to that of _build_graph for test data.
    May be faster than _build_graph for larger ifgs.
    """
    g = graph()
    g.add_nodes(nodes)
    return map(lambda e, w: g.add_edge(edge=e, wt=w), edges, weights)


def mst_matrix_ifgs_only(ifgs):
    '''
    Filter: returns array of independent ifgs from the pixel by pixel MST.

    The MSTs are stripped of connecting edge info, leaving just the ifgs.

    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    '''

    result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix_networkx(ifgs):
        # mst is a dictionary of dates
        if hasattr(mst, 'iteritems'):
            ifg_sub = [ifg_date_lookup(ifgs, d) for d in mst.iteritems()]
            result[(y, x)] = tuple(ifg_sub)
        else:
            result[(y, x)] = mst  # usually NaN
    return result


def mst_matrix_ifg_indices(ifgs):
    '''
    Filter: returns array of independent ifg indices from the pixel by pixel MST.

    The MSTs are stripped of connecting edge info, leaving just the ifgs.
    ifgs: sequence of Ifg objs.

    :param epochs: an EpochList object derived from the ifgs
    '''

    result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix(ifgs):
        if hasattr(mst, 'iteritems'):
            ifg_sub = [ifg_date_index_lookup(ifgs, d) for d in mst.iteritems()]
            result[(y, x)] = tuple(ifg_sub)
        else:
            result[(y, x)] = mst  # usually NaN

    return result


def mst_matrix_as_array(ifgs):
    '''
    Filter: returns array of pixel by pixel MSTs.

    Each pixel contains an MST (with connecting edges etc).

    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    '''

    mst_result = empty(shape=ifgs[0].phase_data.shape, dtype=object)

    for y, x, mst in mst_matrix(ifgs):
        mst_result[y, x] = mst

    return mst_result


# TODO: custom weighting could included with an additional 'weights' arg if some
# other weighting criterion is required later
def mst_matrix(ifgs):
    '''
    Generates/emits MST trees on a pixel-by-pixel basis for the given ifgs.
    :param ifgs: sequence of Ifg objs
    :param epochs: an EpochList object derived from the ifgs
    '''

    # make default MST to optimise result when no Ifg cells in a stack are nans
    epochs = get_epochs(ifgs)
    edges = [(i.master, i.slave) for i in ifgs]
    weights = [i.nan_fraction for i in ifgs]
    g = _build_graph(epochs.dates, edges, weights)
    # TODO: Investigate better algotithms for minimal spanning tree
    dflt_mst = _remove_root_node(minimal_spanning_tree(g))
    # prepare source and dest data arrays
    # TODO: memory efficiencies can be achieved here with tiling
    data_stack = array([i.phase_data for i in ifgs], dtype=float32)

    # create MSTs for each pixel in the ifg data stack
    nifgs = len(ifgs)
    for y, x in product(xrange(i.nrows), xrange(i.ncols)):
        values = data_stack[:, y, x] # vertical stack of ifg values for a pixel
        nancount = sum(isnan(values))

        # optimisations: use pre-created results for all nans/no nans
        if nancount == 0:
            yield y, x, dflt_mst
            continue
        elif nancount == nifgs:
            yield y, x, nan
            continue

        # dynamically modify graph to reuse a single graph: this should avoid
        # repeatedly creating new graph objs & reduce RAM use
        for value, edge, weight in zip(values, edges, weights):
            if not isnan(value):
                if not g.has_edge(edge):
                    g.add_edge(edge, wt=weight)
            else:
                if g.has_edge(edge):
                    g.del_edge(edge)

        mst = _remove_root_node_new(minimal_spanning_tree(g))
        yield y, x, mst


def mst_matrix_networkx(ifgs):
    """
    Generates/emits MST trees on a pixel-by-pixel basis for the given ifgs.
    :param ifgs: sequence of Ifg objs
    """

    # make default MST to optimise result when no Ifg cells in a stack are nans
    epochs = get_epochs(ifgs)
    edges_with_weights_for_networkx = [(i.master, i.slave, i.nan_fraction) for i in ifgs]
    g_nx = _build_graph_networkx(epochs.dates, edges_with_weights_for_networkx)
    T = nx.minimum_spanning_tree(g_nx)

    # TODO: memory efficiencies can be achieved here with tiling
    data_stack = array([i.phase_data for i in ifgs], dtype=float32)

    # create MSTs for each pixel in the ifg data stack
    nifgs = len(ifgs)

    for y, x in product(xrange(i.nrows), xrange(i.ncols)):
        values = data_stack[:, y, x]  # vertical stack of ifg values for a pixel
        nan_count = sum(isnan(values))

        # optimisations: use pre-created results for all nans/no nans
        if nan_count == 0:
            yield y, x, T.edges()
            continue
        elif nan_count == nifgs:
            yield y, x, nan
            continue

        # dynamically modify graph to reuse a single graph: this should avoid
        # repeatedly creating new graph objs & reduce RAM use
        for value, edge in zip(values, edges_with_weights_for_networkx):
            if not isnan(value):
                if not g_nx.has_edge(edge[0], edge[1]):
                    g_nx.add_edge(edge[0], edge[1], weight=edge[2])
            else:
                if g_nx.has_edge(edge[0], edge[1]):
                    g_nx.remove_edges_from(edge)

        yield y, x, T.edges()
