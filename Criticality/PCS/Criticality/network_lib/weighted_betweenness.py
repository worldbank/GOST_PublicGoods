from __future__ import division
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib.pylab import *

from heapq import heappush, heappop
from itertools import count

import os
import pandas as pd
import numpy as np

import networkx as nx
import geopandas as gp

import ema_workbench

from od_prep import od_aggregation

__all__ = ['aon_assignment',
           'probit_assignment',
           'edge_betweenness_centrality',
           'edge_betweenness_subset_od',
           'betweenness_to_df',
           'edge_betweenness_subset_od_ema',
           'ema_betweenness',
           'k_shortest_paths',
           'ksp_edge_betweenness_subset_od',
           'sp_dict_graph_creation',
           'interdiction_single_edge']

def aon_assignment(G, sources, targets, weight, od):
    '''
    Function to do All-or-Nothing assignment on transport network

    Parameters
    ------------
    G: Graph
        Transport network Graph Networkx object that will be analyzed
    sources: list
        List of nodes (integer) that will be used as sources. The integer should correspond to
        node id in G Graph
    targets: list
        List of nodes (integer) that will be used as targets. The integer should correspond to
        node id in G Graph
    weight: str
        String which corresponds to attribute of G Graph's edges that will be used as penalty for each
        edge. In most cases this is defined as 'length' of the edge.
    od: DataFrame
        OD matrix dataframe

    Returns
    ------------
    d: dict
        Dictionary with edge tuple as keys (e.g. (2,3) ) and flow value as values
    '''

    #create empty dict
    d={}

    #iterate over all sources
    for i in range(len(sources)):
        source = sources[i]
        #iterate over all edges
        for j in range(len(targets)):
            target = targets[j]
            #it is assumed that there is no self-loop on the node
            #e.g. there is no flow from node A to node A
            if source != target :
                #determine shortest path between the OD pair
                sp_dijk_all = nx.dijkstra_path(G, source=source, target=target, weight=weight)
                #update the betweenness value of all edges in the shortest path
                flow = od[source][target]
                for j in range(len(sp_dijk_all)-1):
                    lst = [sp_dijk_all[j],sp_dijk_all[j+1]]
                    lst = [min(lst), max(lst)]
                    tup = tuple(lst)
                    if tup in d.keys():
                        d[tup]+=1*flow
                    else:
                        d.update({tup:1*flow})

    #assign 0 to all edges which don't belong to any shortest path
    #at the same time, record all the correct order of edges name
    edges_list = []
    for u,v in G.edges():
        elst = [u,v]
        elst = [min(elst), max(elst)]
        etup = tuple(elst)
        if not etup in d.keys():
            d.update({etup:0})
        tup = tuple([u,v])
        edges_list.append(tup)

    #alter the tuple(u,v) to tuple(v,u) if the order is inconsistent with the original graph's order
    d1 = {}
    for key, val in d.iteritems():
        if not key in edges_list:
            tup = tuple([key[1], key[0]])
            d1.update({tup:val})
        else:
            d1.update({key:val})

    return d1

def probit_assignment(G, sources, targets, weight, od, N=5, sd=10, penalty=0):
    '''
    Function to do stochastic probit assignment on transport network. The weight of the transport network
    is sampled by normal distribution with the original link weight as the mean.

    Parameters
    ------------
    G: Graph
        Transport network Graph Networkx object that will be analyzed
    sources: list
        List of nodes (integer) that will be used as sources. The integer should correspond to
        node id in G Graph
    targets: list
        List of nodes (integer) that will be used as targets. The integer should correspond to
        node id in G Graph
    weight: str
        String which corresponds to attribute of G Graph's edges that will be used as penalty for each
        edge. In most cases this is defined as 'length' of the edge.
    od: DataFrame
        OD matrix dataframe
    N: int
        Number of probit iterations that want to be performed
    sd: int
        Percentage of the link's weight that will be used as standard deviation of the normal distribution (e.g.
        if 10 is inputted, then the standard deviation is 10% of the link's weight). If you don't want to sample
        over the normal distribution, set sd to 0.
    penalty: double
        Penalty that is given to links which have been part of shortest paths set. If set, the value should be higher
        than 1. The intention is to force the code to find distinguished shortest paths between each probit iteration
        by increasing the weight of links that have been part of shortest paths in previous iterations.

    Returns
    ------------
    d: dict
        Dictionary with edge tuple as keys (e.g. (2,3) ) and flow value as values
    '''

    #create empty dict
    d={}

    #create copy of original network to avoid changing the attributes of the original network
    G1 = G.copy()

    #iterate N times
    #in each iteration, sample the link's weight by using normal distribution
    for i in np.arange(N):
        length_dict = {}
        for u,v,data in G1.edges(data=True):
            tup = tuple([u,v])
            if sd > 0:
                length_mean = data[weight]
                stdev = sd/100
                length_sd = length_mean * stdev
                length = np.random.normal(length_mean, length_sd)
                if length < 0:
                    length = 0
            else:
                length = data[weight]
            length_dict.update({tup:length})

        #create a copy of G1 since we want to work the penalty on G1 later
        G2 = G1.copy()

        #set the attribute of G2, we'll work the assignment based on G2's weight information
        nx.set_edge_attributes(G2, weight, length_dict)


        #iterate over all sources
        penalty_list = []
        for i in range(len(sources)):
            source = sources[i]
            #iterate over all edges
            for j in range(len(targets)):
                target = targets[j]
                #it is assumed that there is no self-loop on the node
                #e.g. there is no flow from node A to node A
                if source != target :
                    #determine shortest path between the OD pair
                    sp_dijk_all = nx.dijkstra_path(G2, source=source, target=target, weight=weight)
                    #update the betweenness value of all edges in the shortest path
                    flow = od[source][target]
                    #divide the flow over the number of iteration
                    flow = flow/N
                    for j in range(len(sp_dijk_all)-1):
                        lst = [sp_dijk_all[j],sp_dijk_all[j+1]]
                        lst = [min(lst), max(lst)]
                        tup = tuple(lst)
                        if tup in d.keys():
                            d[tup]+=1*flow
                        else:
                            d.update({tup:1*flow})

                        #if we want to work with penalty, record the shortest paths
                        if penalty > 0:
                            penalty_list.append(tup)
                            tup = tuple([tup[1],tup[0]])
                            penalty_list.append(tup)

        #if work with penalty, update the weight of the links which belong to the shortest paths
        if penalty > 0:
            penalty_dict = {}
            for u,v,data in G1.edges(data=True):
                if tuple([u,v]) in penalty_list:
                    length = data[weight] * penalty
                else:
                    length = data[weight]
                penalty_dict.update({tuple([u,v]):length})

            nx.set_edge_attributes(G1, weight, penalty_dict)

    #assign 0 to all edges which don't belong to any shortest path
    #at the same time, record all the correct order of edges name
    edges_list = []
    for u,v in G.edges():
        elst = [u,v]
        elst = [min(elst), max(elst)]
        etup = tuple(elst)
        if not etup in d.keys():
            d.update({etup:0})
        tup = tuple([u,v])
        edges_list.append(tup)

    #alter the tuple(u,v) to tuple(v,u) if the order is inconsistent with the original graph's order
    d1 = {}
    for key, val in d.iteritems():
        if not key in edges_list:
            tup = tuple([key[1], key[0]])
            d1.update({tup:val})
        else:
            d1.update({key:val})

    return d1

def edge_betweenness_centrality(flow, od):
    '''
    Function to do stochastic probit assignment on transport network. The weight of the transport network
    is sampled by normal distribution with the original link weight as the mean

    Parameters
    ------------
    flow: dict
        Flow dictionary obtained from assignment function (e.g. from aon_assignment or probit_assignment)
    od: DataFrame
        OD matrix dataframe

    Returns
    ------------
    d: dict
        Dictionary with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    '''

    #record the total flow in the network
    totalval = (sum(od.sum()))

    #copy the flow to avoid changing the original flow dictionary
    flow2 = flow.copy()

    #normalize the flow
    for key, val in flow2.items():
        flow2[key] = val / totalval

    return flow2

def edge_betweenness_subset_od(G, sources, targets, weight, od):
    '''
    Old function before betweenness centrality and flow assignment were separated.
    Calculating edge betweenness centrality between only subset of nodes in the network (e.g. between districts)

    Parameters
    ------------
    G: Graph
        Transport network Graph Networkx object that will be analyzed
    sources: list
        List of nodes (integer) that will be used as sources. The integer should correspond to
        node id in G Graph
    targets: list
        List of nodes (integer) that will be used as targets. The integer should correspond to
        node id in G Graph
    weight: str
        String which corresponds to attribute of G Graph's edges that will be used as penalty for each
        edge. In most cases this is defined as 'length' of the edge.
    od: DataFrame
        OD matrix dataframe

    Returns
    ------------
    d: dict
        Dictionary with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    '''

    #create empty dict
    d={}

    #iterate over all sources
    for i in range(len(sources)):
        source = sources[i]
        #iterate over all edges
        for j in range(len(targets)):
            target = targets[j]
            #it is assumed that there is no self-loop on the node
            #e.g. there is no flow from node A to node A
            if source != target :
                #determine shortest path between the OD pair
                sp_dijk_all = nx.dijkstra_path(G, source=source, target=target, weight=weight)
                #update the betweenness value of all edges in the shortest path
                flow = od[source][target]
                for j in range(len(sp_dijk_all)-1):
                    lst = [sp_dijk_all[j],sp_dijk_all[j+1]]
                    lst = [min(lst), max(lst)]
                    tup = tuple(lst)
                    if tup in d.keys():
                        d[tup]+=1*flow
                    else:
                        d.update({tup:1*flow})

    #normalize the betweenness value
    totalval = (sum(od.sum()))
    for key, val in d.items():
        d[key] = val / totalval

    #assign 0 to all edges which don't belong to any shortest path
    for u,v in G.edges():
        elst = [u,v]
        elst = [min(elst), max(elst)]
        etup = tuple(elst)
        if not etup in d.keys():
            d.update({etup:0})

    return d

def betweenness_to_df(gdf,betweenness,betweenness_string):
    '''
    Append betweenness centrality result to the transport network's GeoDataFrame.
    For visualization purpose later.

    Parameters
    ------------
    gdf: GeoDataFrame
        GeoDataFrame (Linestring) of the original transport network
    betweenness: dict
        Dictionary with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    betweenness_string: str
        String of betweenness dictionary's object name

    Returns
    ------------
    gdf_final: GeoDataFrame
        Updated gdf with additional column of betweenness centrality
    betweenness_df: DataFrame
        Betweenness dictionary transformed into dataframe
    '''

    betweenness_df = pd.DataFrame(betweenness.items(), columns=['FromTo_tuple', betweenness_string])

    FromTo_tuple = betweenness_df['FromTo_tuple'].tolist()
    FromTo_tolist = []
    for i in FromTo_tuple:
        odlist = list(i)
        minval = min(odlist)
        maxval = max(odlist)
        val = str(minval) + str(maxval)
        FromTo_tolist.append(val)

    betweenness_df['FromTo'] = FromTo_tolist

    c = []
    for i in range(len(gdf)):
        minval = min([gdf['TNODE_'][i],gdf['FNODE_'][i]])
        maxval = max([gdf['TNODE_'][i],gdf['FNODE_'][i]])
        val = str(minval) + str(maxval)
        c.append(val)
    gdf['FromTo'] = c

    gdf_final = pd.merge(gdf,betweenness_df,on='FromTo',how='outer')

    del gdf_final['FromTo_tuple']

    return gdf_final, betweenness_df

def _shortest_path_record(G, sources, targets, weight):
    '''
    Input:
        G                : Graph Networkx object
        sources, targets : List of nodes sources IDs and nodes targets IDs (e.g. the centroid nodes)
        weight           : Edge data key corresponding to the edge weight
    Output:
        d                : Dict with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    '''

    d={}
    for i in range(len(sources)):
        source = sources[i]
        for j in range(len(targets)):
            target = targets[j]
            if source != target :
                sp_dijk_all = nx.dijkstra_path(G, source=source, target=target, weight=weight)
                od_pair = str(source)+str(target)
                d[od_pair] = (sp_dijk_all, source, target)
    return d

def edge_betweenness_subset_od_ema(G, sp_dict, od):
    d={}
    for key, val in sp_dict.iteritems():
        source = val[1]
        target = val[2]
        sp = val[0]
        flow = od[source][target]
        for j in range(len(sp)-1):
            lst = [sp[j],sp[j+1]]
            lst = [min(lst), max(lst)]
            tup = tuple(lst)
            #the codes below take almost one minute
            if tup in d.keys():
                d[tup]+=1*flow
            else:
                d.update({tup:1*flow})
    totalval = (sum(od.sum()))
    for key, val in d.items():
        d[key] = val / totalval

    for u,v in G.edges():
        elst = [u,v]
        elst = [min(elst), max(elst)]
        etup = tuple(elst)
        if not etup in d.keys():
            d.update({etup:0})
    return d

def ema_betweenness(prod_lists, OD_all_dict, G, sp_dict, **factors_dict):
    OD_final_df = od_aggregation(OD_all_dict, **factors_dict)


    betweenness = edge_betweenness_subset_od_ema(G=G, sp_dict=sp_dict, od=OD_final_df)

    new_d = {}
    for key, val in betweenness.iteritems():
        new_key = str(key[0])+str(key[1])
        new_d[new_key] = val

    return new_d

def k_shortest_paths(G, source, target, k=1, weight='weight'):
    #MAY NOT BE USED ANYMORE
    if source == target:
        return ([0], [[source]])

    length, path = nx.single_source_dijkstra(G, source, target, weight=weight)
    if target not in length:
        raise nx.NetworkXNoPath("node %s not reachable from %s" % (source, target))

    lengths = [length[target]]
    paths = [path[target]]
    c = count()
    B = []
    G_original = G.copy()

    for i in range(1, k):
        for j in range(len(paths[-1]) - 1):
            spur_node = paths[-1][j]
            root_path = paths[-1][:j + 1]

            edges_removed = []
            for c_path in paths:
                if len(c_path) > j and root_path == c_path[:j + 1]:
                    u = c_path[j]
                    v = c_path[j + 1]
                    if G.has_edge(u, v):
                        edge_attr = G.edge[u][v]
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))

            for n in range(len(root_path) - 1):
                node = root_path[n]
                # out-edges
                for u, v, edge_attr in G.edges(node, data=True):
                    G.remove_edge(u, v)
                    edges_removed.append((u, v, edge_attr))

                if G.is_directed():
                    # in-edges
                    for u, v, edge_attr in G.in_edges_iter(node, data=True):
                        G.remove_edge(u, v)
                        edges_removed.append((u, v, edge_attr))

            spur_path_length, spur_path = nx.single_source_dijkstra(G, spur_node, target, weight=weight)
            if target in spur_path and spur_path[target]:
                total_path = root_path[:-1] + spur_path[target]
                total_path_length = _get_path_length(G_original, root_path, weight) + spur_path_length[target]
                heappush(B, (total_path_length, next(c), total_path))

            for e in edges_removed:
                u, v, edge_attr = e
                G.add_edge(u, v, **edge_attr)

        if B:
            (l, _, p) = heappop(B)
            lengths.append(l)
            paths.append(p)
        else:
            break

    return (lengths, paths)

def _get_path_length(G, path, weight='weight'):
    #MAY NOT BE USED ANYMORE
    length = 0
    if len(path) > 1:
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            length += G.edge[u][v].get(weight, 1)

    return length

def _total_cost_sp(G, sources, targets, weight, od):
    '''
    Input:
        G                : Graph Networkx object
        sources, targets : List of nodes sources IDs and nodes targets IDs (e.g. the centroid nodes)
        weight           : Edge data key corresponding to the edge weight
    Output:
        d                : Dict with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    '''
    d={}
    total_cost = 0
    for i in range(len(sources)):
        source = sources[i]
        for j in range(len(targets)):
            target = targets[j]
            if source != target :
                sp_dijk_distance = nx.dijkstra_path_length(G, source=source, target=target, weight=weight)
                flow = od[source][target]
                cost = sp_dijk_distance * flow
                total_cost += cost
                tup=tuple([source,target])
                d.update({tup:cost})

    return total_cost, d

def sp_dict_graph_creation(G, sources, targets, weight):

    sp_dict = _shortest_path_record(G=G, sources = sources, targets = targets, weight=weight)

    edgelist = []
    for edge in list(G.edges()):
        edgelist.append(edge)

    sp_dict_graph = {}
    for key, val in sp_dict.iteritems():
        source = val[1]
        target = val[2]
        tup = tuple([source, target])
        sp_dict_graph.update({tup:[]})
        for j in range(len(val[0])-1):
            test1 = tuple([val[0][j], val[0][j+1]])
            test2 = tuple([val[0][j+1], val[0][j]])
            if test1 in edgelist:
                sp_dict_graph[tup].append(test1)
            if test2 in edgelist:
                sp_dict_graph[tup].append(test2)
    return sp_dict_graph

def interdiction_single_edge(G2, od, weight, sp_dict_graph, sources, targets):
    c = 0
    ff=0
    interdiction_dict = {}
    disconnected_dict = {}
    unsatisfied_demand_dict = {}

    total_cost_base, od_cost_dict = _total_cost_sp(G=G2, sources=sources, targets=targets,
                                                  weight='length', od=od)

    path_in_sp_list = []
    for i in sp_dict_graph.iteritems():
        path_in_sp_list += i[1]

    path_in_sp_list = list(set(path_in_sp_list))

    for i in path_in_sp_list:
        ff += 1
        if ff%200 == 0:
            print(str(ff)+' edges have been interdicted')
        u = i[0]
        v = i[1]
        od_cost_dict2 = od_cost_dict.copy()
        tup = tuple([u,v])
        G = G2.copy()
        G.remove_edge(u,v)
        disconnected = 0
        unsatisfied_demand = 0
        for key, val in sp_dict_graph.iteritems():
            if tup in val:
                try:
                    sp_dijk_distance = nx.dijkstra_path_length(G, source=key[0], target=key[1], weight=weight)
                    flow = od[key[0]][key[1]]
                    cost = sp_dijk_distance * flow
                    od_cost_dict2[key] = cost
                except:
                    sp_dijk_distance = 9999
                    disconnected += 1
                    flow = od[key[0]][key[1]]
                    unsatisfied_demand += flow
        total_cost_new = sum(od_cost_dict2.values())
        cost_increase = (total_cost_new - total_cost_base)/total_cost_base
        unsatisfied_demand = unsatisfied_demand/total_cost_base
        if cost_increase < 0:
            cost_increase = 0
        interdiction_dict.update({tup:cost_increase})
        disconnected_dict.update({tup:disconnected})
        unsatisfied_demand_dict.update({tup:unsatisfied_demand})


    new_interdiction_dict = {}
    for key, val in interdiction_dict.iteritems():
        lst=[key[0], key[1]]
        maxs=max(lst)
        mins=min(lst)
        new_key = str(mins)+str(maxs)
        new_interdiction_dict[new_key] = val

    new_disconnected_dict = {}
    for key, val in disconnected_dict.iteritems():
        lst=[key[0], key[1]]
        maxs=max(lst)
        mins=min(lst)
        new_key = str(mins)+str(maxs)
        new_disconnected_dict[new_key] = val

    new_unsatisfied_demand_dict = {}
    for key, val in unsatisfied_demand_dict.iteritems():
        lst=[key[0], key[1]]
        maxs=max(lst)
        mins=min(lst)
        new_key = str(mins)+str(maxs)
        new_unsatisfied_demand_dict[new_key] = val

    return new_interdiction_dict, new_disconnected_dict, new_unsatisfied_demand_dict

def ksp_edge_betweenness_subset_od(G, sources, targets, weight, od, k):
    '''
    MAY NOT BE USED ANYMORE
    Input:
        G                : Graph Networkx object
        sources, targets : List of nodes sources IDs and nodes targets IDs (e.g. the centroid nodes)
        weight           : Edge data key corresponding to the edge weight
    Output:
        d                : Dict with edge tuple as keys (e.g. (2,3) ) and betweenness value as values
    '''
    d={}
    number=0
    for i in range(len(sources)):
        source = sources[i]
        for j in range(len(targets)):
            target = targets[j]
            if source != target :

                #calculate k-shortest path
                ksp = k_shortest_paths(G = G, source = source, target = target, k = k, weight = weight)

                #store the length of the k-shortest paths
                path_length = ksp[0]
                path_length_set = set(path_length)

                #store total flow between od pair
                flow = od[source][target]

                #calculate logit model for route choice
                #firstly calculate the denominator
                sum_exp = 0
                for i in path_length_set:
                    exp_val = np.exp(-0.05*i)
                    sum_exp += exp_val

                #secondly create list which contains probability of each route
                probability = []
                for i in path_length_set:
                    exp_val = np.exp(-0.05*i)
                    prob = exp_val/sum_exp
                    probability.append(prob)

                #now append the flow*probability to each route
                #iterate for each route
                counter = 0
                for path in path_length_set:
                    index = path_length.index(path)
                    sp = ksp[1][index]
                    for j in range(len(sp)-1):
                        lst = [sp[j],sp[j+1]]
                        lst = [min(lst), max(lst)]
                        tup = tuple(lst)
                        if tup in d.keys():
                            d[tup]+=1*flow*probability[counter]
                        else:
                            d.update({tup:1*flow*probability[counter]})
                    counter += 1

    totalval = (sum(od.sum()))
    for key, val in d.items():
        d[key] = val / totalval

    for u,v in G.edges():
        elst = [u,v]
        elst = [min(elst), max(elst)]
        etup = tuple(elst)
        if not etup in d.keys():
            d.update({etup:0})

    return d