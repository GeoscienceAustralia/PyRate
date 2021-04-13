from typing import List
from pyrate.core.logger import pyratelogger as log


def dfs(graph, marked, n, vert, start, count, loop, all_loops):
    """
    Python Program to count cycles of length n in a given graph. Number of vertices V
    This is an adaptation of https://www.geeksforgeeks.org/print-all-the-cycles-in-an-undirected-graph/
    """
    V = graph.shape[0]

    # mark the vertex vert as visited
    marked[vert] = True

    # if the path of length (n-1) is found
    if n == 0:

        # mark vert as un-visited to make it usable again.
        marked[vert] = False

        # Check if vertex vert can end with vertex start
        if graph[vert][start] == 1:
            count += 1
            all_loops.append(loop)
            return count
        else:
            return count

    # For searching every possible path of length (n-1)
    for i in range(V):
        if (not marked[i]) and (graph[vert][i] == 1):
            next_path = loop[:]
            next_path.append(i)
            # DFS for searching path by decreasing length by 1
            count = dfs(graph, marked, n - 1, i, start, count, next_path, all_loops)

    # marking vert as unvisited to make it usable again.
    marked[vert] = False
    return count


def find_cycles(graph, loop_length):
    """Counts cycles of length N in an undirected and connected graph"""

    count, all_loops = find_all_cycled_of_length(graph, loop_length)
    deduped_loops = __dedupe_loops(all_loops)
    return deduped_loops


def find_all_cycled_of_length(graph, loop_length):
    log.info(f"finding loops of length {loop_length}")
    V = graph.shape[0]
    all_loops = []
    # all vertex are marked un-visited initially.
    marked = [False] * V
    # Searching for cycle by using v-n+1 vertices
    count = 0
    for i in range(V - (loop_length - 1)):
        count = dfs(graph, marked, loop_length - 1, i, i, count, [i], all_loops)

        # ith vertex is marked as visited and will not be visited again.
        marked[i] = True
    log.debug(f"Total possible number of loops of length {loop_length} is {count}")
    return count, all_loops


def __dedupe_loops(simple_cycles: List[List]) -> List:
    """
    Deduplication of loops with same members.

    For example: for nodes 1, 2, 3 in the graph below

        [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]

    Loops of length 3 are 0->1-> 2 and 0->2->1. We only retain 1 of these loops after dedupe.

    Example 2:
    For the following graph:
            [
                [0, 1, 1, 1],
                [1, 0, 1, 1],
                [1, 1, 0, 1],
                [1, 1, 1, 0],
            ]

    Pictorially this is the network
    1 ---------2
    | *      * |
    |    *     |
    |  *    *  |
    4----------3

    Loops of length 4 are:
        [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1]]

    After deduplication we will only retain [0, 1, 2, 3]

    """
    seen_sc_sets = set()
    filtered_sc = []
    for sc in simple_cycles:
        loop = sc[:]
        sc.sort()
        sc = tuple(sc)
        if sc not in seen_sc_sets:
            seen_sc_sets.add(sc)
            filtered_sc.append(loop)
    log.info(f"Selected loops after deduplication is {len(filtered_sc)}")
    return filtered_sc
