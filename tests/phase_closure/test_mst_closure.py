from datetime import date
from pathlib import Path
import numpy as np
import pytest
from pyrate.constants import PYRATEPATH
from pyrate.core.phase_closure.mst_closure import (
    find_closed_loops, Edge, SignedWeightedEdge, SignedEdge, setup_edges,
    add_signs_and_weights_to_loops, sort_loops_based_on_weights_and_date, WeightedLoop,
    find_signed_closed_loops
)

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs


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
    assert isinstance(signed_loops[0], WeightedLoop)
    swe = signed_loops[0].loop[0]
    assert isinstance(swe, SignedWeightedEdge)
    assert isinstance(swe.signed_edge, SignedEdge)
    assert isinstance(swe.edge, Edge)
    assert isinstance(swe.first, date)
    assert isinstance(swe.second, date)


def test_sort_loops_based_on_weights_and_date(signed_loops):
    weighted_loops = sort_loops_based_on_weights_and_date(signed_loops)
    assert len(weighted_loops) == 541
    # order
    weights = [w.weight for w in weighted_loops]
    earliest_dates = [w.earliest_date for w in weighted_loops]
    assert np.all(np.diff(weights) >= 0)

    for i, (w, d, wl) in enumerate(zip(weights[:-1], earliest_dates[:-1], weighted_loops[:-1])):
        sub_list = [d]
        if w == weights[i+1]:
            sub_list.append(earliest_dates[i+1])
        if len(sub_list) > 1:
            tds = np.array([td.days for td in np.diff(sub_list)])
            assert np.all(tds >= 0)  # assert all dates are increasing for same weights


def test_add_signs_and_weights_to_loops(geotiffs):
    """also tests find_signed_closed_loops"""
    all_edges = setup_edges(geotiffs)
    all_loops = find_closed_loops(all_edges)
    loops1 = add_signs_and_weights_to_loops(all_loops, all_edges)

    all_edges = setup_edges(geotiffs)
    all_loops = find_closed_loops(all_edges)
    loops2 = add_signs_and_weights_to_loops(all_loops, all_edges)

    compare_loops(loops1, loops2)


def compare_loops(loops1, loops2):
    for i, (l1, l2) in enumerate(zip(loops1, loops2)):
        assert l1.weight == l2.weight
        for ll1, ll2 in zip(l1.loop, l2.loop):
            assert ll1.weight == ll2.weight
            assert ll1.first == ll2.first
            assert ll1.second == ll2.second
            assert ll1.sign == ll2.sign
    assert i == len(loops1) - 1
    assert i == 540


def test_find_signed_closed_loops(geotiffs):
    loops1 = find_signed_closed_loops(geotiffs)
    loops2 = find_signed_closed_loops(geotiffs)
    compare_loops(loops1, loops2)


def test_sort_loops_based_on_weights_and_date_2(geotiffs):
    loops1 = find_signed_closed_loops(geotiffs)
    loops2 = find_signed_closed_loops(geotiffs)
    sorted_loops1 = sort_loops_based_on_weights_and_date(loops1)
    sorted_loops2 = sort_loops_based_on_weights_and_date(loops2)
    compare_loops(sorted_loops1, sorted_loops2)
