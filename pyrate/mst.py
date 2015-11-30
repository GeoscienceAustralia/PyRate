"""
Minimum Spanning Tree functionality for PyRate.

Contains functions to calculate MST using interferograms.

.. codeauthors:: Ben Davies, Sudipta Basak
"""

from itertools import product
from numpy import array, nan, isnan, float32, empty

from pyrate.algorithm import get_all_epochs, master_slave_ids
from pyrate.algorithm import ifg_date_lookup, get_epochs
from pyrate.algorithm import ifg_date_index_lookup

import networkx as nx

# TODO: may need to implement memory saving row-by-row access
# TODO: document weighting by either Nan fraction OR variance


# TODO: this does not appear to be used anywhere
def default_mst(ifgs):
    '''
    Returns default MST dict for the given Ifgs. The MST is calculated using a
    weighting based on the number of incoherent cells in the phase band.

    :param noroot: True removes the PyGraph default root node from the result.
    '''

    edges_with_weights_for_networkx = [(i.master, i.slave, i.nan_fraction)
                                       for i in ifgs]
    g_nx = _build_graph_networkx(edges_with_weights_for_networkx)
    mst = nx.minimum_spanning_tree(g_nx).edges()
    return mst


def _build_graph_networkx(edges_with_weights):
    """
    Convenience graph builder function: returns a new graph obj.
    """
    g = nx.Graph()
    g.add_weighted_edges_from(edges_with_weights)
    return g


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

    for y, x, mst in mst_matrix_networkx(ifgs):
        if isinstance(mst, list):
            ifg_sub = [ifg_date_index_lookup(ifgs, d) for d in mst]
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
    edges_with_weights_for_networkx = [(i.master, i.slave, i.nan_fraction)
                                       for i in ifgs]
    g_nx = _build_graph_networkx(edges_with_weights_for_networkx)
    T = nx.minimum_spanning_tree(g_nx)
    edges = T.edges()

    # TODO: memory efficiencies can be achieved here with tiling
    data_stack = array([i.phase_data for i in ifgs], dtype=float32)

    # create MSTs for each pixel in the ifg data stack
    nifgs = len(ifgs)

    for y, x in product(xrange(i.nrows), xrange(i.ncols)):
        values = data_stack[:, y, x]  # vertical stack of ifg values for a pixel
        nan_count = sum(isnan(values))

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
        for value, edge in zip(values, edges_with_weights_for_networkx):
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
