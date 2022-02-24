#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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

from typing import List
from pyrate.core.logger import pyratelogger as log


def dfs(graph, marked, n, vert, start, count, loop, all_loops):
    """
    Depth First Search algorithm to count loops of length n in a given graph.
    Adapted from https://www.geeksforgeeks.org/print-all-the-cycles-in-an-undirected-graph/
    """
    # Number of vertices
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


def find_loops(graph, loop_length):
    """Counts and collects loops of a defined length in an undirected and connected graph"""
    V = graph.shape[0]
    all_loops = []
    # all vertex are marked un-visited initially.
    marked = [False] * V
    # Searching for loop by using v-n+1 vertices
    count = 0
    for i in range(V - (loop_length - 1)):
        count = dfs(graph, marked, loop_length - 1, i, i, count, [i], all_loops)

        # ith vertex is marked as visited and will not be visited again.
        marked[i] = True
    log.debug(f"Total possible number of loops of length {loop_length} is {count}")
    return count, all_loops


def dedupe_loops(loops: List[List]) -> List:
    """
    Deduplication of loops with same members.

    For example: for nodes 1, 2, 3 in the graph below

        [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]

    Loops of length 3 are 0->1->2 and 0->2->1. We only retain 1 of these loops after dedupe.

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
    seen_sets = set()
    filtered = []
    for l in loops:
        loop = l[:]
        l.sort()
        l = tuple(l)
        if l not in seen_sets:
            seen_sets.add(l)
            filtered.append(loop)
    return filtered
