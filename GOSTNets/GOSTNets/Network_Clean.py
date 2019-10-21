#the cleaning network part
import os, sys, time
sys.path.append(r'/home/vagrant/repos/GOST_PublicGoods/GOSTNets/GOSTNets')
import GOSTnet as gn
import importlib
import networkx as nx
import osmnx as ox
from shapely.ops import unary_union
from shapely.wkt import loads
from shapely.geometry import LineString, MultiLineString, Point


def CleanNetwork(G, wpath = '', country='', UTM={'init': 'epsg:3857'}, WGS = {'init': 'epsg:4326'}, junctdist = 50, verbose = False):
    ''' Topologically simplifies an input graph object by collapsing junctions and removing interstital nodes
    REQUIRED - G: a graph object containing nodes and edges. edges should have a property 
                  called 'Wkt' containing geometry objects describing the roads
                wpath: the write path - a drive directory for inputs and output
                country: this parameter allows for the sequential processing of multiple countries
                UTM: the epsg code of the projection, in metres, to apply the junctdist
    OPTIONAL - junctdist: distance within which to collapse neighboring nodes. simplifies junctions. 
                Set to 0.1 if not simplification desired. 50m good for national (primary / secondary) networks
                verbose: if True, saves down intermediate stages for dissection
    '''
    # Squeezes clusters of nodes down to a single node if they are within the snapping tolerance
    a = gn.simplify_junctions(G, UTM, WGS, junctdist)

    # ensures all streets are two-way
    a = gn.add_missing_reflected_edges(a)
    
    #save progress
    if verbose is True: 
        gn.save(a, 'a', wpath)
    
    # Finds and deletes interstital nodes based on node degree
    b = gn.custom_simplify(a)
    
    # rectify geometry
    for u, v, data in b.edges(data = True):
        if type(data['Wkt']) == list:
                data['Wkt'] = gn.unbundle_geometry(data['Wkt'])
    
    # save progress
    if verbose is True: 
        gn.save(b, 'b', wpath)
    
    # For some reason CustomSimplify doesn't return a MultiDiGraph. Fix that here
    c = gn.convert_to_MultiDiGraph(b)

    # This is the most controversial function - removes duplicated edges. This takes care of two-lane but separate highways, BUT
    # destroys internal loops within roads. Can be run with or without this line
    c = gn.remove_duplicate_edges(c)

    # Run this again after removing duplicated edges
    c = gn.custom_simplify(c)

    # Ensure all remaining edges are duplicated (two-way streets)
    c = gn.add_missing_reflected_edges(c)
    
    # save final
    if verbose:
        gn.save(c, '%s_processed' % country, wpath)
    
    print('Edge reduction: %s to %s (%d percent)' % (G.number_of_edges(), 
                                               c.number_of_edges(), 
                                               ((G.number_of_edges() - c.number_of_edges())/G.number_of_edges()*100)))
    return c

def largestSubgraph(G_copy):
    ''' Extract the largest sub graph from an input networkx object
    
    INPUT
    G_copy [networkx digraph]
    
    RETRUNS
    [networkx digraph]
    '''
    list_of_Gs = list((nx.strongly_connected_component_subgraphs(G_copy)))
    list_length = list(len(i) for i in list_of_Gs)
    m = max(list_length)
    t = [i for i, j in enumerate(list_length) if j == m][0]
    max_G = list_of_Gs[t]
    return(max_G)
