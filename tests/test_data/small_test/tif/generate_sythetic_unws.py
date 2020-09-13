import networkx as nx
from pyrate.core.mst import _minimum_spanning_edges_from_mst

g = nx.Graph()
g.add_weighted_edges_from(edges_with_weights)
return g

