from __future__ import division
import pandas as pd
import numpy as np

import networkx as nx
import geopandas as gp
from osmnx.utils import make_str

from shapely.geometry import LineString, Point

from osmnx_simplify_overwrite import simplify_graph

from weighted_betweenness import probit_assignment

def prepare_centroids_network2(centroid, network):
    '''
    Take transport network and centroids shapefiles as inputs
    then returns a geodataframe of the transport network with
    indicated centroid nodes

    Parameters
    ------------
    centroid: str
        string of centroid shapefile's address+filename
    network: str
        string of network shapefile's address+name

    Returns
    ------------
    gdf_points: GeoDataFrame
        geodataframe (Points) of centroids shapefile
    gdf_node_pos: GeoDataFrame
        geodataframe (Points) of nodes obtained from all links in the network shapefile.
        The IsCentroid information is attached if a node is the closest node from a centroid.
    gdf: GeoDataFrame
        geodataframe (LineString) of the original network, containing information about the start
        node (FNODE) and end node (TNODE). The IsCentroid information is attached as well.
    '''

    #read the centroid shapefile into geodataframe
    gdf_points = gp.read_file(centroid)

    #read the network shapefile into geodataframe
    gdf = prepare_gdf_network(network)

    #take all nodes from the network geodataframe into dataframe
    df_node_pos = gdf[['Start_pos', 'FNODE_', 'road']].rename(columns={'Start_pos': 'pos', 'FNODE_': 'Node', 'road': 'road' }).append(
                gdf[['End_pos', 'TNODE_', 'road']].rename(columns={'End_pos': 'pos', 'TNODE_': 'Node', 'road': 'road' }))

    #drop all duplicate nodes
    df_node_pos = df_node_pos.drop_duplicates(subset='Node')

    #change the column name
    df_node_pos.columns = ['geometry', 'Node', 'road']

    #add column of POINT type for the geometry
    df_node_pos['pointgeo'] = [Point(xy) for xy in df_node_pos.geometry]

    #reindex the dataframe
    df_node_pos.index = range(len(df_node_pos))

    #save the longitude (x) and latitude(y) separately
    xy = np.array(df_node_pos['geometry'].tolist())
    x = [xy[i,0] for i in range(len(xy))]
    y = [xy[i,1] for i in range(len(xy))]
    df_node_pos['x'] = x
    df_node_pos['y'] = y

    #create geodataframe of the network points from dataframe
    gdf_node_pos = gp.GeoDataFrame(df_node_pos, crs=gdf.crs, geometry=df_node_pos.pointgeo)
    gdf_node_pos['osmid'] = gdf_node_pos.index

    #reference the Node ID of the network to the centroids by selecting the nearest node from the centroid points
    gdf_points['Node'] = gdf_points.geometry.apply(lambda g: gdf_node_pos.iloc[gdf_node_pos.distance(g).idxmin()].Node)
    OD = gdf_points['Node'].tolist()
    gdf_node_pos['IsCentroid'] = gdf_node_pos.Node.apply(lambda g: 1 if g in OD else 0)

    #adding Centroid information to the gdf
    gdf['IsCentroid1'] = gdf.TNODE_.apply(lambda g: 1 if g in OD else 0)
    gdf['IsCentroid2'] = gdf.FNODE_.apply(lambda g: 1 if (g in OD) else 0)
    gdf['IsCentroid'] = gdf['IsCentroid1'] + gdf['IsCentroid2']
    del gdf['IsCentroid1']
    del gdf['IsCentroid2']

    #create unique osmid for the network LineString GeoDataFrame
    gdf['osmid'] = gdf.index.map(lambda x: x + 10000)

    return gdf_points, gdf_node_pos, gdf


__all__ = ['prepare_gdf_network',
           'prepare_centroids_network',
           'gdf_to_simplified_multidigraph',
           'multigraph_to_graph',
           'graph_to_df',
           'prepare_adm_background',
           'create_link_capacity']

def prepare_gdf_network(network):
    '''
    Converting transport network shapefile into GeoDataFrame

    Parameters
    ------------
    network: str
        string of network shapefile's address+filename

    Returns
    ------------
    gdf: GeoDataFrame
        geodataframe of network with linestring, coordinate of start position, and
        coordinate of end position recorded
    '''

    # Load network shapefile into GeoDataFrame
    gdf = gp.read_file(network)

    # !!! Add column capacity for min max cut flow algorithm
    # gdf['capacity'] = gdf['RD_CLASS']

    # shapefile needs to include minimal: geometry linestring and the length computed (e.g. in QGIS)
    if 'length' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if not gdf.geometry.map(lambda x: type(x) ==  LineString).all():
        s_invalid_geo = gdf.geometry[gdf.geometry.map(lambda x: type(x) ==  LineString)]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))

    # Compute the start- and end-position based on linestring
    gdf['Start_pos'] = gdf.geometry.apply(lambda x: x.coords[0])
    gdf['End_pos'] = gdf.geometry.apply(lambda x: x.coords[-1])

    # Create Series of unique nodes and their associated position
    s_points = gdf.Start_pos.append(gdf.End_pos).reset_index(drop=True)
    s_points = s_points.drop_duplicates()

    # Add index of start and end node of linestring to geopandas DataFrame
    df_points = pd.DataFrame(s_points, columns=['Start_pos'])
    df_points['FNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='Start_pos', how='inner')

    df_points = pd.DataFrame(s_points, columns=['End_pos'])
    df_points['TNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='End_pos', how='inner')

    return gdf

def prepare_centroids_network(centroid, network):
    '''
    Take transport network and centroids shapefiles as inputs
    then returns a geodataframe of the transport network with
    indicated centroid nodes

    Parameters
    ------------
    centroid: str
        string of centroid shapefile's address+filename
    network: str
        string of network shapefile's address+name

    Returns
    ------------
    gdf_points: GeoDataFrame
        geodataframe (Points) of centroids shapefile
    gdf_node_pos: GeoDataFrame
        geodataframe (Points) of nodes obtained from all links in the network shapefile.
        The IsCentroid information is attached if a node is the closest node from a centroid.
    gdf: GeoDataFrame
        geodataframe (LineString) of the original network, containing information about the start
        node (FNODE) and end node (TNODE). The IsCentroid information is attached as well.
    '''

    #read the centroid shapefile into geodataframe
    gdf_points = gp.read_file(centroid)

    #read the network shapefile into geodataframe
    gdf = prepare_gdf_network(network)

    #take all nodes from the network geodataframe into dataframe
    df_node_pos = gdf[['Start_pos', 'FNODE_']].rename(columns={'Start_pos': 'pos', 'FNODE_': 'Node' }).append(
                gdf[['End_pos', 'TNODE_']].rename(columns={'End_pos': 'pos', 'TNODE_': 'Node' }))

    #drop all duplicate nodes
    df_node_pos = df_node_pos.drop_duplicates(subset='Node')

    #change the column name
    df_node_pos.columns = ['geometry', 'Node']

    #add column of POINT type for the geometry
    df_node_pos['pointgeo'] = [Point(xy) for xy in df_node_pos.geometry]

    #reindex the dataframe
    df_node_pos.index = range(len(df_node_pos))

    #save the longitude (x) and latitude(y) separately
    xy = np.array(df_node_pos['geometry'].tolist())
    x = [xy[i,0] for i in range(len(xy))]
    y = [xy[i,1] for i in range(len(xy))]
    df_node_pos['x'] = x
    df_node_pos['y'] = y

    #create geodataframe of the network points from dataframe
    gdf_node_pos = gp.GeoDataFrame(df_node_pos, crs=gdf.crs, geometry=df_node_pos.pointgeo)
    gdf_node_pos['osmid'] = gdf_node_pos.index

    #reference the Node ID of the network to the centroids by selecting the nearest node from the centroid points
    gdf_points['Node'] = gdf_points.geometry.apply(lambda g: gdf_node_pos.iloc[gdf_node_pos.distance(g).idxmin()].Node)
    OD = gdf_points['Node'].tolist()
    gdf_node_pos['IsCentroid'] = gdf_node_pos.Node.apply(lambda g: 1 if g in OD else 0)

    #adding Centroid information to the gdf
    gdf['IsCentroid1'] = gdf.TNODE_.apply(lambda g: 1 if g in OD else 0)
    gdf['IsCentroid2'] = gdf.FNODE_.apply(lambda g: 1 if (g in OD) else 0)
    gdf['IsCentroid'] = gdf['IsCentroid1'] + gdf['IsCentroid2']
    del gdf['IsCentroid1']
    del gdf['IsCentroid2']

    #create unique osmid for the network LineString GeoDataFrame
    gdf['osmid'] = gdf.index.map(lambda x: x + 10000)

    return gdf_points, gdf_node_pos, gdf

def gdf_to_simplified_multidigraph(gdf_node_pos, gdf, undirected = True, simplify=True, sums=['length']):
    '''
    Simplifying transport network (in GeoDataFrame format) by removing all nodes which are
    neither intersections, end/start nodes, nor centroids. This reduces the computation time
    needed to conduct any analysis later.

    Parameters
    ------------
    gdf_node_pos: GeoDataFrame
        geodataframe (Points) of nodes obtained from all links in the network shapefile.
        The IsCentroid information is attached if a node is the closest node from a centroid.
    gdf: GeoDataFrame
        geodataframe (LineString) of the original transport network.

    Returns
    ------------
    G2: MultiGraph, MultiDiGraph
        Multi(Di)Graph Networkx object of simplified transport network. Multi(Di)Graph network type
        is chosen because the graph simplification algorithm only works with this type of network.
    '''

    #create a MultiDiGraph object
    G2 = nx.MultiDiGraph(crs=gdf.crs)

    #create nodes on the MultiDiGraph
    for index, row in gdf_node_pos.iterrows():
        c =  {'x': row.x, 'y': row.y, 'IsCentroid' : row.IsCentroid, 'ID' : row.Node, 'osmid': row.osmid}
        G2.add_node(row.Node, **c)

    #create bidirectional edges on top of the MultiDiGraph nodes
    #based on the FNODE and TNODE information of the transport network GeoDataFrame
    for index, row in gdf.iterrows():
        dict_row  = row.to_dict()
        if 'geometry' in dict_row: del dict_row['geometry']
        G2.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)

    #swap the FNODE and the TNODE since we want to make bidirectional graph
    gdf.rename(columns={'Start_pos': 'End_pos',
                       'End_pos': 'Start_pos',
                       'FNODE_': 'TNODE_',
                       'TNODE_': 'FNODE_', }, inplace=True)

    #do the iteration again
    for index, row in gdf.iterrows():
        dict_row  = row.to_dict()
        if 'geometry' in dict_row: del dict_row['geometry']
        G2.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)

    #simplify the MultiDiGraph using OSMNX's overwritten function
    if simplify:
        G2 = simplify_graph(G_=G2, sums=sums)

    #make a name
    G2.graph['name'] = 'graph'

    if undirected:
        G2 = G2.to_undirected()
    return G2

def multigraph_to_graph(G):
    '''
    Change Multi(Di)Graph object to Graph. Graph is undirected, simple graph type without parallel edges
    (while Multi(Di)Graph may have parallel edges). This code removes duplicate edges.

    Parameters
    ------------
    G: MultiGraph, MultiDiGraph
        Multi(Di)Graph Networkx object of simplified transport network

    Returns
    ------------
    G2_new_tograph: Graph
        Graph Networkx object
    '''

    #create arbitrary Graph object
    G2_new_tograph = nx.Graph()

    #create dummy Graph as a mean to indicate duplicated edges
    G_dummy = nx.Graph()

    #transform the nodes into Graph, preserving all attributes
    for u,v in G.nodes(data=True):
        G2_new_tograph.add_node(u, **v)
        G_dummy.add_node(u, **v)

    #transform the edges into Graph, preserving all attributes
    c = []
    for u,v,data in G.edges(data=True):
        d = (u,v)
        #check if the edge that connects (u,v) exists in the graph
        if not d in c:
            G2_new_tograph.add_edge(u,v,**data)
            G_dummy.add_edge(u,v,**data)
            c.append(d)
        #else replace the old edge with the new edge if the new edge has longer length
        else:
            for edge in G_dummy.edges(data=True):
                e = (edge[0], edge[1])
                if e == d:
                    if data['length'] > edge[2]['length']:
                        G2_new_tograph.remove_edge(u,v)
                        G2_new_tograph.add_edge(u,v,**data)
                        G_dummy.remove_edge(u,v)
                        G_dummy.add_edge(u,v,**data)

    return G2_new_tograph

def graph_to_df(G2_simplified):
    '''
    Change Graph Networkx object back to GeoDataFrame. This helps for visualization purpose,
    as GeoDataFrame has more flexibility in displaying the transport network.

    Parameters
    ------------
    G2_simplified: Graph
        (Simplified) Graph Networkx object

    Returns
    ------------
    gdf_edges: GeoDataFrame
        GeoDataFrame (Linestring) of the Graph Networkx object
    '''

    #get undirected Graph from MultiDiGraph
    G2_simplified3 = G2_simplified.copy()

    #create a copy for safer operation
    G_save = G2_simplified3.copy()

    #create dictionaries of nodes from the undirected Graph
    nodes = {node:data for node, data in G_save.nodes(data=True)}

    #create GeoDataFrame of nodes
    gdf_nodes = gp.GeoDataFrame(nodes).T

    #change the CRS (coordinate reference system) into EPSG:4326
    gdf_nodes.crs = {'proj': 'longlat', 'ellps': 'WGS84', 'datum': 'WGS84'}

    #define the geometry attribute of the GeoDataFrame
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    gdf_nodes = gdf_nodes.drop(['x', 'y'], axis=1)

    #ensure the osmid is in integer
    gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(np.int64)

    #remove all nodes that do not have geometry information
    for col in [c for c in gdf_nodes.columns if not c == 'geometry']:
        gdf_nodes[col] = gdf_nodes[col].fillna('').map(make_str)

    #create list of edges
    edges = []
    for u, v, data in G_save.edges(data=True):
        edge_details = data
        edge_details['FNODE_'] = u
        edge_details['TNODE_'] = v

        # if edge doesn't already have a geometry attribute, create one now
        if not 'geometry' in data:
            point_u = Point((G_save.node[u]['x'], G_save.node[u]['y']))
            point_v = Point((G_save.node[v]['x'], G_save.node[v]['y']))
            edge_details['geometry'] = LineString([point_u, point_v])

        edges.append(edge_details)

    #create GeoDataFrame of edges
    gdf_edges = gp.GeoDataFrame(edges)

    #change the CRS into EPSG:4326
    gdf_edges.crs = {'proj': 'longlat', 'ellps': 'WGS84', 'datum': 'WGS84'}

    #if any element of the GeoDataFrame contains more than one information (due to the graph simplification)
    for i in gdf_edges.columns:
        #select only one element which is most relevant
        gdf_edges[i] = gdf_edges[i].apply(lambda x: _last_list_element(x))

    #delete all irrelevant columns
    del gdf_edges['End_pos']
#     del gdf_edges['FNODE_']
#     del gdf_edges['TNODE_']
#     del gdf_edges['osmid']
    del gdf_edges['Start_pos']
    # del gdf_edges['capacity']

    return gdf_edges

def _last_list_element(x):
    #utility for multidigraph_to_shapefile function
    #take the last element of a list
    if type(x) == list:
        x = x[-1]
    return x

def prepare_adm_background(adm_csv, adm_shp, csv_column_list):
    '''
    Preparing geospatial administrative area background in GeoDataFrame. Merge various socioeconomic
    information from another csv. The csv and the shp should have at least one column with identical values
    to map to each other.

    Parameters
    ------------
    adm_csv: str
        string of socioeconomic csv file address+filename
    adm_shp: str
        string of administrative area csv file address+filename
    csv_column_list: list
        list of string of column names from the adm_csv file that want to be added into the
        resulting GeoDataFrame

    Returns
    ------------
    district_gdf2: GeoDataFrame
        GeoDataFrame (Polygon) of administrative area and its corresponding socioeconomic data
    '''

    #read district data statistics
    district_df = pd.read_csv(adm_csv)

    #read adm level 2 district shapefile
    district_gdf = gp.read_file(adm_shp)

    #extract only the intended columns
    district_df2 = district_df[csv_column_list]

    #change the 'code' column into 'HASC_2' so that it can be merged with the shp file
    district_df2.rename(columns={'Code':'HASC_2'}, inplace=True)

    #combine gdf and df
    district_gdf2 = pd.merge(district_gdf, district_df2, on='HASC_2')

    return district_gdf2

def create_link_capacity(G, item1, item2='length', calctype='multiplication'):
    '''
    Preparing capacity of a link for unimodal transport network (i.e. road or waterway or railway separately).
    This function (currently) only perform simple multiplication or division between the two items
    in order to generate capacity attribute.

    Parameters
    ------------
    G: Graph
        Networkx Graph object of a unimodal transport network
    item1: str
        string of the Graph's data attribute that want to be used as the first component of capacity calculation
    item2: str
        string of the Graph's data attribute that want to be used as the secoond component of capacity calculation

    Returns
    ------------
    G1: Graph
        Networkx Graph object of a unimodal transport network with capacity attribute embedded in its edges
    '''

    capacity_dict = {}
    G1 = G.copy()

    for u,v,data in G.edges(data=True):
        if type(data[item1]) != list:
            component1 = data[item1]
        else:
            component1 = min(data[item1])

        if type(data[item2]) != list:
            component2 = data[item2]
        else:
            component2 = min(data[item2])

        edge = tuple([u,v])
        if calctype == 'multiplication':
            capacity = component1 * component2
        elif calctype == 'division':
            capacity = component1 / component2

        capacity_dict.update({edge:capacity})

    nx.set_edge_attributes(G1, 'lcapacity', capacity_dict)

    return G1