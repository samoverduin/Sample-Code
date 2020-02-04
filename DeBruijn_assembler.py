#!/usr/bin/env python

"""
Author: Sam Overduin
Student nr: 920208637010
Script to: find eulerian paths in spectrum of l-mers
"""
# Import statements
import random


def count_inout_degree(graph_dct):
    """Returns a dct of vertices as values and (in_degree,out_degree) as values

    :param graph_dct: directed graph as key with lst of keys
    :return: dct of vertice keys and (in_degree, out_degree) as values.
    """
    out_dct = {key:0 for key in graph_dct.keys()}
    in_dct = {}
    for out_node, in_lst in graph_dct.items():
        out_dct[out_node] = len(in_lst)
        for in_node in in_lst:
            try:
                in_dct[in_node] += 1
            except:
                in_dct[in_node] = 1

    inout_dct = {}
    in_keys = list(in_dct.keys())
    out_keys = list(out_dct.keys())
    nodes = set(list(in_keys+out_keys))
    for node in nodes:
        # every node represented in in_dct and out_dct
        if node not in in_keys:
            in_dct[node] = 0
        if node not in out_keys:
            out_dct[node] = 0
        
        inout_dct[node] = (in_dct[node], out_dct[node])
    return inout_dct


def is_eulerian(graph_dct):
    """Returns a bool for if graph is eulerian

    :param graph_dct: directed graph as key with lst of keys
    :return: True/False if all vertices are balanced
    """
    inout_dct = count_inout_degree(graph_dct)
    for key, tup in inout_dct.items():
        in_cnt, out_cnt = tup
        if in_cnt != out_cnt:
            return False
    return True


def has_eulerian_path(graph_dct, return_start_end_nodes=False):
    """Returns a bool for if graph has a eulerian path

    :param graph_dct: directed graph as key with lst of keys
    :param return_start_end_nodes: bool to return start+end nodes
    :return: True/False if max 2 vertices have |in_deg-out_deg| = 1,
        node start (vertice with 1 less in than out),
        node end (vertice with 1 more in than out)

    """
    deg_diffs = 0
    inout_dct = count_inout_degree(graph_dct)
    start_node = None
    end_node = None
    has_path = True

    for key, tup in inout_dct.items():
        in_cnt, out_cnt = tup
        if in_cnt - out_cnt == 1:
            if not end_node:
                end_node = key
            else:
                has_path, end_node = (False, None)
                break
        elif in_cnt - out_cnt == -1:
            if not start_node:
                start_node = key
            else:
                has_path, start_node = (False, None)
                break
        elif abs(in_cnt - out_cnt) > 1:
             has_path, start_node, end_node = (False, None, None)

    if return_start_end_nodes:
        return has_path, start_node, end_node
    else:
        return has_path


def make_euler_graph(spect):
    """ Creates a (l-1)-mer graph where the edges correspond to every l-mer.

    :param spect: a lst of l-mer sequences (spectrum)
    :return: directed graph as key with lst of keys
    """
    assert(len(set([len(x) for x in spect])) == 1),\
        "spect elements must be same length"

    graph_dct = {}
    for seq in spect:
        try:
            graph_dct[seq[:-1]].append(seq[1:])
        except KeyError:
            graph_dct[seq[:-1]] = [seq[1:]]
    return graph_dct


def find_eulerian_cycle(graph_dct, start_point=None):
    """ Returns a eulerian cycle containing all edges of the graph

    :param graph_dct: directed graph as key with lst of keys
    :param start_point: optional key in graph_dct to start from
    :return: list of vertices in order.
    """
    edges_left = 0
    for key, to_lst in graph_dct.items():
        edges_left += len(to_lst)

    if not start_point:
        start_point = random.choice(list(graph_dct.keys()))

    path = []
    left_paths = {}
    cur_point = start_point
    new_point = None

    while edges_left != 0 and new_point != start_point:
        # randomized edge to follow
        new_point = graph_dct[cur_point].pop(
            random.randint(0, len(graph_dct[cur_point])-1))
        path.append((cur_point, new_point))
        left_paths[cur_point] = len(graph_dct[cur_point])
        cur_point = new_point
        edges_left -= 1
    # Check unused edges
    for key, v in left_paths.items():
        if v > 0:
            extra_path = find_eulerian_cycle(graph_dct, start_point=key)
            # find where to cut open and place extra path in between
            i = [start for start, end in path].index(key)
            path = path[:i] + extra_path + path[i:]
    return path


def find_eulerian_path(graph):
    """ Returns a eulerian path as lst of keys from graph_dct
        by visiting all edges once

    :param graph: directed graph as key with lst of keys
    :return: lst of keys of eulerian path.
    """
    # copy due to aliasing and popping of stuff in graph
    graph_dct = {k:v[:] for k, v in graph.items()}

    if not is_eulerian(graph_dct):
        has_path, start, end = has_eulerian_path(graph_dct,
                                                 return_start_end_nodes=True)
        if not has_path:
            print("graph_dct does not contain Eulerian path.")
            return None
        # make euler cycle:
        try:
            graph_dct[end].append(start)
        except KeyError:
            graph_dct[end] = [start]

        path_lst = find_eulerian_cycle(graph_dct, start_point=start)
        # remove last one because not actually a cycle
        path_lst = path_lst[:-1]
    else:
        path_lst = find_eulerian_cycle(graph_dct)

    assert len(path_lst) == len(set(path_lst)), "duplicate edges in path"

    # Take all first nodes from edges + add final end node (2nd line)
    node_lst = [start for start, end in path_lst]
    node_lst.append(path_lst[-1][1])

    return node_lst


def main():
    # Eulerian graph
    graph_822 = {'J': ['D'], 'D': ['C'], 'I': ['H'], 'H': ['F'], 'F': ['G', 'E'],
                 'C': ['I', 'A'], 'G': ['J'], 'E': ['A'], 'A': ['F', 'B'], 'B': ['C']}
    # 3-mers created from sequence (not in order)
    s = ['ATG', 'TGG', 'TGC', 'GTG', 'GGC', 'GCA', 'GCG', 'CGT']

    print('Eulerian path through 822:')
    print(find_eulerian_path(graph_822))

    print("\nCreate graph from k-mers s & find path through the graph")
    s_graph = make_euler_graph(s)
    seq_path = find_eulerian_path(s_graph)
    print('path:', seq_path)
    tot_seq = seq_path[0]
    for i in range(1, len(seq_path)):
        tot_seq += seq_path[i][-1:]
    print('sequence:', tot_seq)


if __name__ == "__main__":
    main()
