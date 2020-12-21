from pathlib import Path
import pytest
from pyrate.constants import PYRATEPATH
from pyrate.core.phase_closure.mst_closure import (
    find_closed_loops, Edge, SignedWeightedEdge, SignedEdge, setup_edges,
    add_signs_and_weights_to_loops,
)

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    return [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]


@pytest.fixture
def all_loops(geotiffs):
    edges = setup_edges(geotiffs)
    loops = find_closed_loops(edges)
    assert len(loops) == 541
    return loops


@pytest.fixture
def edges(geotiffs):
    all_edges = setup_edges(geotiffs)
    return all_edges


@pytest.fixture
def signed_loops(all_loops, edges):
    loops = add_signs_and_weights_to_loops(all_loops, edges)
    return loops


@pytest.fixture(params=[True, False])
def weight(request):
    return request.param


def test_setup_edges(geotiffs):
    edges = setup_edges(geotiffs)
    assert len(edges) == len(geotiffs) == 30
    assert isinstance(edges[0], Edge)


def test_associate_ifgs_with_loops(signed_loops, geotiffs):
    assert len(geotiffs) == 30
    assert len(signed_loops) == 541
    assert isinstance(signed_loops[0][0], SignedWeightedEdge)
    assert isinstance(signed_loops[0][0][0], SignedEdge)
    assert isinstance(signed_loops[0][0][0][0], Edge)
