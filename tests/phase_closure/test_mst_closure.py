from pathlib import Path
import pytest
from pyrate.constants import PYRATEPATH
from pyrate.core.phase_closure.mst_closure import (
    find_closed_loops, Edge, WeightedEdge, setup_edges,
    filter_loops_to_increasing_sequence_loops,
    associate_ifgs_with_loops
)

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    return [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]


@pytest.fixture
def all_loops(geotiffs):
    weighted_edges = setup_edges(geotiffs, weighted=True)
    loops = find_closed_loops(weighted_edges)
    assert len(loops) == 5468  # TODO: improve
    return loops


@pytest.fixture
def loops(all_loops):
    increaings_loops = filter_loops_to_increasing_sequence_loops(all_loops)
    assert len(increaings_loops) == 75
    return increaings_loops


def test_setup_edges():
    pass


def test_associate_ifgs_with_loops(loops, geotiffs):
    edges = setup_edges(geotiffs, weighted=False)
    loops_to_ifgs = associate_ifgs_with_loops(loops, edges)
    assert len(loops_to_ifgs) == 75
