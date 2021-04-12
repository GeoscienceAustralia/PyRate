from typing import List


def dfs(graph, marked, n, vert, start, count, path, all_loops):
    """
    Python Program to count cycles of length n in a given graph. Number of vertices V
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
            count = count + 1
            all_loops.append(path)
            return count
        else:
            return count

    # For searching every possible path of length (n-1)
    for i in range(V):
        if (not marked[i]) and (graph[vert][i] == 1):
            next_path = path[:]
            next_path.append(i)
            # DFS for searching path by decreasing length by 1
            count = dfs(graph, marked, n - 1, i, start, count, next_path, all_loops)

    # marking vert as unvisited to make it usable again.
    marked[vert] = False
    return count


def find_cycles(graph, n):
    """Counts cycles of length N in an undirected and connected graph"""
    V = graph.shape[0]
    all_loops = []

    # all vertex are marked un-visited initially.
    marked = [False] * V

    # Searching for cycle by using v-n+1 vertices
    count = 0
    for i in range(V - (n - 1)):
        count = dfs(graph, marked, n - 1, i, i, count, [i], all_loops)

        # ith vertex is marked as visited and will not be visited again.
        marked[i] = True
    print('count of loops: ', len(all_loops), count)
    deduped_loops = __dedupe_loops(all_loops)
    return count/2, deduped_loops


def __dedupe_loops(simple_cycles: List[List]) -> List:
        """ also discard loops when the loop members are the same"""
        seen_sc_sets = set()
        filtered_sc = []
        for sc in simple_cycles:
            loop = sc[:]
            sc.sort()
            sc = tuple(sc)
            print('tuple: ', sc)
            print(seen_sc_sets)
            if sc not in seen_sc_sets:
                seen_sc_sets.add(sc)
                print(sc, 'adding:', loop)
                filtered_sc.append(loop)
            else:
                print(sc, 'not adding:', loop)

        return filtered_sc
