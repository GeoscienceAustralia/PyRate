import numpy as np
from pyrate.core.phase_closure.collect_loops import find_cycles


def test_collect_loops():
    graph = np.array(
        [
            [0, 1, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [0, 1, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [0, 1, 0, 1, 0]
        ]
    )

    n = 4
    count, all_loops = find_cycles(graph, n)
    assert count == 3
    np.testing.assert_array_equal(all_loops, [[0, 1, 2, 3], [0, 1, 4, 3], [1, 2, 3, 4]])
