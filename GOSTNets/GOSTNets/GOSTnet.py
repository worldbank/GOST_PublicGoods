import peartree as pt
import peartree.graph as ptg
print('peartree version: %s ' % pt.__version__)
import networkx as nx
print('networkx version: %s ' % nx.__version__)
import matplotlib as mpl
print('matplotlib version: %s ' % mpl.__version__)
import osmnx as ox
print('osmnx version: %s ' % ox.__version__)
import os, sys
import pandas as pd, geopandas as gpd
from functools import partial
import pyproj
from shapely.ops import transform, linemerge
from shapely.wkt import loads
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import linemerge, unary_union
import time
import numpy as np
from collections import Counter
import warnings

def node_gdf_from_graph(G, crs = {'init' :'epsg:4326'}, attr_list = None):

    #### Function for generating GeoDataFrame from Graph ####
    # REQUIRED: a graph object G
    # OPTIONAL: crs - projection of format {'init' :'epsg:4326'}. Defaults to WGS84.
    #           note: defining crs of the data, does NOT reproject to this crs
    #           attr_list: list of the keys which you want to be moved over to the GeoDataFrame, if not all
    # RETURNS: a geodataframe of the node objects in the graph
    # -------------------------------------------------------------------------#
    nodes = []
    keys = []

    if attr_list is None:
        for u, data in G.nodes(data = True):
            keys.append(list(data.keys()))
        flatten = lambda l: [item for sublist in l for item in sublist]
        keys = list(set(flatten(keys)))
        attr_list = keys

    if 'geometry' in attr_list:
        non_geom_attr_list = attr_list
        non_geom_attr_list.remove('geometry')
    else:
        non_geom_attr_list = attr_list

    z = 0

    for u, data in G.nodes(data=True):

        if 'geometry' not in attr_list and 'x' in attr_list and 'y' in attr_list :
            try:
                new_column_info = {
                'node_ID': u,
                'geometry': Point(data['x'], data['y']),
                'x': data['x'],
                'y': data['y']}
            except:
                print((u, data))
        else:
            try:
                new_column_info = {
                'node_ID': u,
                'geometry': data['geometry'],
                'x':data['geometry'].x,
                'y':data['geometry'].y}
            except:
                print((u, data))
        for i in non_geom_attr_list:
            try:
                new_column_info[i] = data[i]
            except:
                pass

        nodes.append(new_column_info)
        z += 1

    nodes_df = pd.DataFrame(nodes)
    nodes_df = nodes_df[['node_ID',*non_geom_attr_list,'geometry']]
    nodes_df = nodes_df.drop_duplicates(subset=['node_ID'], keep='first')
    nodes_gdf = gpd.GeoDataFrame(nodes_df, geometry=nodes_df.geometry, crs = crs)

    return nodes_gdf

def edge_gdf_from_graph(G, crs = {'init' :'epsg:4326'}, attr_list = None, geom_col = 'geometry'):

    #### Function for generating GeoDataFrame from Graph ####
    # REQUIRED: a graph object G
    # OPTIONAL: crs - projection of format {'init' :'epsg:4326'}. Defaults to WGS84
    #           attr_list: list of the keys which you want to be moved over to the GeoDataFrame
    # RETURNS: a GeoDataFrame object of the edges in the graph
    # -------------------------------------------------------------------------#


    edges = []
    keys = []

    if attr_list is None:
        for u, v, data in G.edges(data = True):
            keys.append(list(data.keys()))
        flatten = lambda l: [item for sublist in l for item in sublist]
        keys = list(set(flatten(keys)))
        if geom_col in keys:
            keys.remove(geom_col)
        attr_list = keys

    for u, v, data in G.edges(data=True):

        if geom_col in data:
            # if it has a geometry attribute (a list of line segments), add them
            # to the list of lines to plot
            geom = data[geom_col]

        else:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[u]['x']
            y1 = G.nodes[u]['y']
            x2 = G.nodes[v]['x']
            y2 = G.nodes[v]['y']
            geom = LineString([(x1, y1), (x2, y2)])

        new_column_info = {
            'stnode':u,
            'endnode':v,
            'geometry':geom}

        for i in attr_list:
            try:
                new_column_info[i] = data[i]
            except:
                pass

        edges.append(new_column_info)

    edges_df = pd.DataFrame(edges)
    edges_df = edges_df[['stnode','endnode',*attr_list,'geometry']]
    edges_gdf = gpd.GeoDataFrame(edges_df, geometry = edges_df.geometry, crs = crs)

    return edges_gdf

def snap_points_to_graph(G, points, response = None, geomcol = 'geometry', connection_threshold = 5000):

    #### Function for generating GeoDataFrame from Graph ####
    # REQUIRED: a GeoDataFrame of point objects (points_gdf)
    #           a Graph object or geodataframe
    # OPTIONAL: response: return result in different formats - dataframe, list, or
    #           list of unique nearest nodes
    #           geomcol: specify a different column name to the default for the
    #           geometry column of the input GeoDataframe. Useful for GeoDataFrames
    #           that include multiple columns with shapely geometry info
    #           connection threshold: stop considering nodes further than the connection_threshold. Default 5km.
    # RETURNS:  an augmented version of the input GeoDataFrame with the node_ID of
    #           the nearest nodes to the points in the graph
    # Note:     ensure any GeoDataFrames are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#

    node_df_G1 = node_gdf_from_graph(G)

    if type(points) != gpd.geodataframe.GeoDataFrame:
        raise ValueError('points variable must be of type GeoDataFrame!')

    nn,nl = [], []

    for i, row in points.iterrows():

        pointobj = points[geomcol].loc[i]

        point = (pointobj.x, pointobj.y)

        nearest_nodes = get_nearest_nodes(node_df_G1, point, connection_threshold = connection_threshold)

        try:
            nrst_node = nearest_nodes.end_node.loc[nearest_nodes.length.idxmin()]
            nrst_node_dist = nearest_nodes.length.loc[nearest_nodes.length.idxmin()]
        except:
            nrst_node = None
            nrst_node_dist = None

        nn.append(nrst_node)
        nl.append(nrst_node_dist)

    points['Nearest_node_ID'] = nn
    points['Nearest_node_dist'] = nl

    if response is None:
        return points
    elif response == 'list':
        return list(points['Nearest_node'])
    elif response == 'unique_list':
        return list(set(nn))
    else:
        ValueError('response parameter not recongized!')

def graph_nodes_intersecting_polygon(G, polygons, crs = None):

    #### Function for generating GeoDataFrame from Graph ####
    # REQUIRED: a GeoDataFrame containing one or more polygons
    #           a Graph object or geodataframe
    # RETURNS:  a list of the nodes intersecting the polygons
    # Note:     ensure any GeoDataFrames are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#



    if type(G) == nx.classes.multidigraph.MultiDiGraph:
        graph_gdf = node_gdf_from_graph(G)

    elif type(G) == gpd.geodataframe.GeoDataFrame:
        graph_gdf = G
    else:
        raise ValueError('Expecting a graph or geodataframe for G!')

    if type(polygons) != gpd.geodataframe.GeoDataFrame:
        raise ValueError('Expecting a geodataframe for polygon(s)!')

    if crs != None and graph_gdf.crs != crs:
            graph_gdf = graph_gdf.to_crs(crs)

    if crs != None and polygons.crs != crs:
            polygons = polygons.to_crs(crs)

    if polygons.crs != graph_gdf.crs:
        raise ValueError('crs mismatch detected! aborting process')

    aggs = []
    for poly in polygons.geometry:

        def chck(x, poly):
            if poly.contains(x):
                return 1
            else:
                return 0

        graph_gdf['intersecting'] = graph_gdf['geometry'].apply(lambda x: chck(x, poly))
        aggs.append(list(graph_gdf['node_ID'].loc[graph_gdf['intersecting'] == 1]))

    aggs = [j for i in aggs for j in i]
    aggs = list(set(aggs))
    return aggs

def graph_edges_intersecting_polygon(G, polygons, mode = 'contains', crs = None):

    #### Function for identifying intersecting edges of a graph with polygon(s) ####
    # REQUIRED: a GeoDataFrame containing one or more polygons
    #           a Graph object
    #           mode - a string, either 'contains' or 'intersecting'
    # RETURNS:  a list of the edges intersecting the polygons
    # Note:     ensure any GeoDataFrames are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#

    if type(G) == nx.classes.multidigraph.MultiDiGraph:
        node_graph_gdf = node_gdf_from_graph(G)
        edge_graph_gdf = edge_gdf_from_graph(G)
    else:
        raise ValueError('Expecting a graph or geodataframe for G!')

    if type(polygons) != gpd.geodataframe.GeoDataFrame:
        raise ValueError('Expecting a geodataframe for polygon(s)!')

    if crs != None and node_graph_gdf.crs != crs:
            node_graph_gdf = node_graph_gdf.to_crs(crs)

    if crs != None and polygons.crs != crs:
            polygons = polygons.to_crs(crs)

    if polygons.crs != node_graph_gdf.crs:
        raise ValueError('crs mismatch detected! aborting process')

    intersecting_nodes = graph_nodes_intersecting_polygon(node_graph_gdf, polygons, crs)

    if mode == 'contains':
        edge_graph_gdf = edge_graph_gdf.loc[(edge_graph_gdf.stnode.isin(intersecting_nodes)) &
                                 (edge_graph_gdf.endnode.isin(intersecting_nodes))]
    elif mode == 'intersects':
        edge_graph_gdf = edge_graph_gdf.loc[(edge_graph_gdf.stnode.isin(intersecting_nodes)) |
                                 (edge_graph_gdf.endnode.isin(intersecting_nodes))]

    return edge_graph_gdf

def sample_raster(G, tif_path, property_name = 'RasterValue'):

    #### Function for attaching raster values to corresponding graph nodes ####
    # REQUIRED: a graph containing one or more nodes
    #           a raster or path to a tif
    #           a property name for the value of the raster attached to the node
    # RETURNS:  a graph
    # Note:     ensure any GeoDataFrames / graphs are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#

    import rasterio

    if type(G) == nx.classes.multidigraph.MultiDiGraph or type(G) == nx.classes.digraph.DiGraph:
        pass
    else:
        raise ValueError('Expecting a graph or geodataframe for G!')

    # generate dictionary of {node ID: point} pairs
    try:
        list_of_nodes = {}
        for u, data in G.nodes(data=True):
            list_of_nodes.update({u:(data['x'], data['y'])})
    except:
        raise ValueError('loading point geometry went wrong. Ensure data dict includes x, y values!')

    # load raster
    try:
        dataset = rasterio.open(os.path.join(tif_path))
    except:
        raise ValueError('Expecting a path to a .tif file!')

    # create list of values
    raster_values = list(dataset.sample(list_of_nodes.values()))
    raster_values = [x[0] for x in raster_values]

    # generate new dictionary of {node ID: raster values}
    ref = dict(zip(list_of_nodes.keys(), raster_values))

    # load new values onto node data dictionary
    for u, data in G.nodes(data=True):
        data[property_name] = ref[u]

    return G

def generate_isochrones(G, origins, thresh, weight = None, stacking = False):

    #### Function for generating isochrones from one or more graph nodes ####
    # REQUIRED: G - a graph containing one or more nodes
    #           origins - a list of node IDs that the isochrones are to be generated from
    #           thresh - the time threshold for the calculation of the isochrone
    # OPTIONAL: weight - name of edge weighting for calculating 'distances'. For isochrones, should be
    #           time expressed in seconds. Defaults to time expressed in seconds.
    #           stacking - if True, returns number of origins that can be reached from that node. If false, max = 1
    # RETURNS:  the original graph with a new data property for the nodes and edges included in the isochrone
    # Note:     ensure any GeoDataFrames / graphs are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#

    if type(origins) == list and len(origins) >= 1:
        pass
    else:
        raise ValueError('Ensure isochrone centers (origins object) is a list containing at least one node ID!')

    ddict = list(G.nodes(data = True))[:1][0][1]

    if weight == None:
        if 'time' not in ddict.keys():
            raise ValueError('need "time" key in edge value dictionary!')
        else:
            weight = 'time'

    sub_graphs = []
    for node in origins:
        sub_graphs.append(nx.ego_graph(G, node, thresh, distance = weight))

    reachable_nodes = []
    for graph in sub_graphs:
        reachable_nodes.append(list(graph.nodes))

    reachable_nodes = [j for i in reachable_nodes for j in i]

    if stacking == False:

        reachable_nodes = set(reachable_nodes)

        for u, data in G.nodes(data=True):
            if u in reachable_nodes:
                data[thresh] = 1
            else:
                data[thresh] = 0

    elif stacking == True:

        reachable_nodes = Counter(reachable_nodes)

        for u, data in G.nodes(data=True):
            if u in reachable_nodes:
                data[thresh] = reachable_nodes[u]
            else:
                data[thresh] = 0
    else:
        raise ValueError('stacking must either be True or False!')

    return G

def make_iso_polys(G, origins, trip_times, edge_buff=25, node_buff=50, infill=False, weight = None, crs = None):

    default_crs = {'init':'epsg:4326'}

    ddict = list(G.nodes(data = True))[:1][0][1]

    if type(origins) == list and len(origins) >= 1:
        pass
    else:
        raise ValueError('Ensure isochrone centers ("origins" object) is a list containing at least one node ID!')

    if weight == None:
        if 'time' not in ddict.keys():
            raise ValueError('need "time" key in edge value dictionary!')
        else:
            weight = 'time'

    isochrone_polys, nodez, tt = [], [], []

    for trip_time in sorted(trip_times, reverse=True):

        for _node_ in origins:

            subgraph = nx.ego_graph(G, _node_, radius = trip_time, distance = weight)

            node_points = [Point((data['x'], data['y'])) for node, data in subgraph.nodes(data=True)]
            nodes_gdf = gpd.GeoDataFrame({'id': subgraph.nodes()}, geometry=node_points, crs = default_crs)
            nodes_gdf = nodes_gdf.set_index('id')

            edge_lines = []
            for n_fr, n_to in subgraph.edges():
                f = nodes_gdf.loc[n_fr].geometry
                t = nodes_gdf.loc[n_to].geometry
                edge_lines.append(LineString([f,t]))

            edge_gdf = gpd.GeoDataFrame({'geoms':edge_lines}, geometry = 'geoms', crs = default_crs)

            if crs != None and nodes_gdf.crs != crs:
                nodes_gdf = nodes_gdf.to_crs(crs)
                edge_gdf = edge_gdf.to_crs(crs)

            n = nodes_gdf.buffer(node_buff).geometry
            e = edge_gdf.buffer(edge_buff).geometry

            all_gs = list(n) + list(e)

            new_iso = gpd.GeoSeries(all_gs).unary_union

            # If desired, try and "fill in" surrounded
            # areas so that shapes will appear solid and blocks
            # won't have white space inside of them

            if infill:
                new_iso = Polygon(new_iso.exterior)

            isochrone_polys.append(new_iso)
            tt.append(trip_time)
            nodez.append(str(_node_))

    gdf = gpd.GeoDataFrame({'geometry':isochrone_polys,'thresh':tt,'nodez':_node_}, crs = crs, geometry = 'geometry')

    return gdf

def convert_network_to_time(G, distance_tag, graph_type = 'drive', road_col = 'highway', speed_dict = None, walk_speed = 4.5, factor = 1):

    #### Function for adding a time value to edge dictionaries ####
    # REQUIRED: G - a graph containing one or more nodes
    #           distance_tag - the key in the dictionary for the field currently containing a distance in meters
    # OPTIONAL: road_col - key for the road type in the edge data dictionary
    #           graph_type - flags network type
    #           speed_dict - speed dictionary to use. If not supplied, reverts to defaults
    #           walk_speed - specify a walkspeed in km/h
    #           factor - allows you to scale up / down distances if saved in
    #                    a unit other than metres
    # RETURNS:  the original graph with a new data property for the edges called 'time
    # Note:     ensure any GeoDataFrames / graphs are in the same projection
    #           before using function, or pass a crs
    # -------------------------------------------------------------------------#

    ## TODO ##
    # deal with graphs with multiple edges between node pairs

    if type(G) == nx.classes.multidigraph.MultiDiGraph or type(G) == nx.classes.digraph.DiGraph:
        pass
    else:
        raise ValueError('Expecting a graph or geodataframe for G!')

    G_adj = G.copy()

    for u, v, data in G_adj.edges(data=True):

        orig_len = data[distance_tag] * factor

        # Note that this is a MultiDiGraph so there could
        # be multiple indices here, I naively assume this is not
        # the case
        data['length'] = orig_len

        # get appropriate speed limit
        if graph_type == 'walk':
            speed = walk_speed

        elif graph_type == 'drive':

            if speed_dict == None:
                speed_dict = {
                'residential': 20,  # kmph
                'primary': 40, # kmph
                'primary_link':35,
                'motorway':50,
                'motorway_link': 45,
                'trunk': 40,
                'trunk_link':35,
                'secondary': 30, # kmph
                'secondary_link':25,
                'tertiary':30,
                'tertiary_link': 25,
                'unclassified':20
                }
            highwayclass = data[road_col]

            if type(highwayclass) == list:
                highwayclass = highwayclass[0]

            if highwayclass in speed_dict.keys():
                speed = speed_dict[highwayclass]
            else:
                speed = speed_dict['residential']

        else:
            raise ValueError('Expecting either a graph_type of "walk" or "drive"!')

        # perform conversion
        kmph = (orig_len / 1000) / speed
        in_seconds = kmph * 60 * 60
        data['time'] = in_seconds

        # And state the mode, too
        data['mode'] = graph_type

    return G_adj

def bind_graphs(G1,G2,name, exempt_nodes, connection_threshold = 50, speed = 4.5, verbose = True):

    # Terminate this process early if either graph is empty
    if (G1.number_of_nodes() == 0) or (G2.number_of_nodes() == 0) :
        return pd.DataFrame({'stop_id': [],
                             'to_nodes': [],
                             'edge_costs': []})

    # First, we need a DataFrame representation of the nodes in the graph
    node_df_G1 = node_gdf_from_graph(G1)
    node_df_G2 = node_gdf_from_graph(G2)

    # Remove all nodes that are part of the new additions to the graph
    if len(exempt_nodes) > 0:
        node_df_G2 = node_df_G2[~node_df_G2.index.isin(exempt_nodes)]

    nn = []

    i, j = 1, 0
    ten_pct = int(len(node_df_G2) / 10)
    for i, row in node_df_G2.iterrows():

        sid = str(row.node_ID)
        full_sid = ptg.nameify_stop_id(name, sid)

        # Ensure that each value is typed correctly prior to being
        # fed into the nearest node method
        lon = float(row.x)
        lat = float(row.y)
        point = (lon, lat)

        nearest_nodes = get_nearest_nodes(node_df_G1,
                                          point,
                                          connection_threshold,
                                          exempt_id=full_sid)

        # Iterate through series results and add to output
        nearest_nodes['start_node'] = sid
        nearest_nodes['start_point'] = Point(point)
        nearest_nodes['mode'] = 'network_binding'
        nearest_nodes = nearest_nodes.loc[nearest_nodes.length < connection_threshold]

        nn.append(nearest_nodes)

        if i % ten_pct == 0 and verbose == True:
            print('    finished binding %d percent of nodes' % (j*10))
            j+= 1
        i += 1

    nearest_nodes = pd.concat(nn)

    nearest_nodes['end_point'] = nearest_nodes.apply(lambda x: Point(x.end_point_x, x.end_point_y), axis = 1)
    nearest_nodes = nearest_nodes[['start_node','start_point','end_node','end_point','length']]
    nearest_nodes['geometry'] = nearest_nodes.apply(lambda x: LineString((x.start_point, x.end_point)), axis = 1)

    Gnew = G1.copy()
    G2_nodes = list(G2.nodes(data = True))
    Gnew.add_nodes_from(G2_nodes)
    G2_edges = list(G2.edges(data = True))
    Gnew.add_edges_from(G2_edges)

    edge_list = []

    for i, row in nearest_nodes.iterrows():

        orig_len = row.length

        kmph = (orig_len / 1000) / speed
        in_seconds = kmph * 60 * 60

        e = (row.start_node, row.end_node, {'length':orig_len, 'time':in_seconds})
        f = (row.end_node, row.start_node, {'length':orig_len, 'time':in_seconds})

        edge_list.append(e)
        edge_list.append(f)

    Gnew.add_edges_from(edge_list)

    return Gnew

def great_circle_vec(lat1: float,
                     lng1: float,
                     lat2: float,
                     lng2: float,
                     earth_radius: float=6371009.0) -> float:
    """
    Vectorized function to calculate the great-circle distance between two
    points or between vectors of points.
    Please note that this method is copied from OSMnx method of the same name,
    which can be accessed here:
    https://github.com/gboeing/osmnx/blob/
    b32f8d333c6965a0d2f27c1f3224a29de2f08d55/osmnx/utils.py#L262
    Parameters
    ----------
    lat1 : float or array of float
    lng1 : float or array of float
    lat2 : float or array of float
    lng2 : float or array of float
    earth_radius : numeric
        radius of earth in units in which distance will be returned (default is
        meters)
    Returns
    -------
    distance : float
        distance or vector of distances from (lat1, lng1) to (lat2, lng2) in
        units of earth_radius
    """
    import warnings

    phi1 = np.deg2rad(90 - lat1)
    phi2 = np.deg2rad(90 - lat2)

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)

    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) \
           + np.cos(phi1) * np.cos(phi2))

    # Ignore warnings during this calculation because numpy warns it cannot
    # calculate arccos for self-loops since u==v
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        arc = np.arccos(cos)

    # Return distance in units of earth_radius
    distance = arc * earth_radius
    return distance

def get_nearest_nodes(df_orig: pd.DataFrame,
                      point,
                      connection_threshold: float,
                      exempt_id: str=None):
    # This method breaks out a portion of a similar method from
    # OSMnx's get_nearest_node; source:
    #   https://github.com/gboeing/osmnx/blob/
    #   b32f8d333c6965a0d2f27c1f3224a29de2f08d55/osmnx/utils.py#L326

    # Make a copy of the DataFrame to prevent mutation outside of function
    df = df_orig.copy()

    if exempt_id is not None:
        df.index = df.index.astype(str)
        mask = ~(df.index == exempt_id)
        df = df[mask]

    # Add second column of reference points
    df['reference_y'] = point[1]
    df['reference_x'] = point[0]

    # TODO: OSMnx supports euclidean as well, for now we have a stumped
    #       version of this same function

    # Ensure each vectorized series is typed correctly
    ref_ys = df['reference_y'].astype(float)
    ref_xs = df['reference_x'].astype(float)
    ys = df['y'].astype(float)
    xs = df['x'].astype(float)

    # Calculate distance vector using great circle distances (ie, for
    # spherical lat-long geometries)
    distances = great_circle_vec(lat1=ref_ys,
                                 lng1=ref_xs,
                                 lat2=ys,
                                 lng2=xs)

    # Filter out nodes outside connection threshold

    nearest_nodes = pd.DataFrame({'length':distances,
                                  'end_node':df.node_ID,
                                 'end_point_x':df.x,
                                  'end_point_y':df.y})

    # Return filtered series
    return nearest_nodes

def example_edge(G, n):
    i = list(G.edges(data = True))[:n]
    for j in i:
        print(j)

def example_node(G, n):
    i = list(G.nodes(data = True))[:n]
    for j in i:
        print(j)

def calculate_OD(G, origins, destinations, fail_value):
    #### Function for generating an origin: destination matrix  ####
    # REQUIRED: G - a graph containing one or more nodes
    #           fail_value - the value to return if the trip cannot be completed (implies some sort of disruption / disconnected nodes)
    #           origins - a list of the node IDs to treat as origins points
    #           destinations - a list of the node IDs to treat as destination points
    # RETURNS:  a numpy matrix of format OD[o][d] = shortest time possible
    # -------------------------------------------------------------------------#

    OD = np.zeros((len(origins), len(destinations)))

    for o in range(0, len(origins)):
        origin = origins[o]

        for d in range(0, len(destinations)):
            destination = destinations[d]

            try:
                shortest_time = nx.shortest_path_length(G, source=origin, target=destination, weight='time')
                OD[o][d] = shortest_time
            except:
                OD[o][d] = fail_value

    return OD

def disrupt_network(G, property, thresh, fail_value):
    #### Function for disrupting a graph ####
    # REQUIRED: G - a graph containing one or more nodes and one or more edges
    #           property - the element in the data dictionary for the edges to test
    #           thresh - values of data[property] above this value are disrupted
    #           fail_value - the data['time'] property is set to this value to simulate the removal of the edge
    # RETURNS:  a modified graph with the edited 'time' attribute
    # -------------------------------------------------------------------------#

    G_copy = G.copy()

    broken_nodes = []

    for u, data in G_copy.nodes(data = True):

        if data[property] > thresh:

            broken_nodes.append(u)
    print(broken_nodes)

    for u, v, data in G.edges(data = True):

        if u in broken_nodes or v in broken_nodes:

            data['time'] = fail_value

    return G_copy

def randomly_disrupt_network(G, edge_frac, fail_value):

    #### Function for randomly disurpting a network ####
    # REQUIRED: G - a graph containing one or more nodes and one or more edges
    #           edge_frac - the percentage of edges to destroy. Interger rather than decimal
    #           fail_value - the data['time'] property is set to this value to simulate the removal of the edge
    # RETURNS:  a modified graph with the edited 'time' attribute
    #           the list of edge IDs randomly chosen for destruction
    # NOTE:     requires the graph to have an 'edge_id' value in the edge data dictionary. This DOES NOT have to be unique.
    # -------------------------------------------------------------------------#

    edgeid = []

    for u,v, data in G.edges(data = True):
        edgeid.append(data['edge_id'])

    num_to_destroy = math.floor(len(edgeid) / 2 * (edge_frac / 100))

    destroy_list = list(np.random.randint(low = 0, high = max(edgeid), size = [num_to_destroy]))

    G_adj = G.copy()

    for u, v, data in G_adj.edges(data = True):
        if data['edge_id'] in destroy_list:
            data['time'] = fail_value

    return G_adj, destroy_list

def test_connectivity(G, target_node):
    #### Function for testing which nodes are currently connected to a network ####
    # REQUIRED: G - a graph containing one or more nodes and one or more edges
    #           edge_frac - the percentage of edges to destroy. Interger rather than decimal
    #           fail_value - the data['time'] property is set to this value to simulate the removal of the edge
    # RETURNS:  a modified graph with the edited 'time' attribute
    #           the list of edge IDs randomly chosen for destruction
    # NOTE:     requires the graph to have an 'edge_id' value in the edge data dictionary. This DOES NOT have to be unique.
    # -------------------------------------------------------------------------#

    answers, successnodes = [], []
    worked, failed = 0, 0

    for u in G.nodes:
        v = target_node
        try:
            calc = nx.shortest_path_length(G, source=u, target=v, weight='time')
            answers.append({'origin':u,'destination':v,'spath':calc})
            successnodes.append(u)
            worked +=1
        except:
            answers.append({'origin':u,'destination':v,'spath':None})
            failed +=1

    #### Add attribute for whether node passed connectivity trial
    for u, data in G.nodes(data = True):
        if u in successnodes:
            data['test'] = 'ok'
        else:
            data['test'] = 'fail'

    print('connected nodes: %s' % worked)
    print('broken nodes: %s' % failed)

def gravity_demand(G, origins, destinations, weight, maxtrips = 100, dist_decay = 1, fail_value = 99999999999):
    #### Function for generating a gravity-model based demand matrix ####
    # REQUIRED: G - a graph containing one or more nodes and one or more edges
    #           origins - a list of node IDs
    #           destinations - a list of node IDs
    #           weight - the gravity weighting of the nodes in the model, e.g. population
    # OPTIONAL: fail_value - the data['time'] property is set to this value to simulate the removal of the edge
    #           maxtrips - normalize the number of trips in the resultant function to this number of trip_times
    #           dist_decay - parameter controlling the aggresion of discounting based on distance
    # RETURNS:  a numpy array describing the demand between o and d in terms of number of trips
    # NOTE:     1 trip will always be returned between an origins and a destination, even if weighting would otherewise be 0
    # -------------------------------------------------------------------------#

    maxtrips = 100
    dist_decay = 1

    demand = np.zeros((len(origins), len(destinations)))

    shortest_time = Calculate_OD(G, origins, destinations, fail_value)

    for o in range(0, len(origins)):
        for d in range(0, len(destinations)):
            if origins == destinations and o == d:
                demand[o][d] = 0
            else:
                normalized_dist = shortest_time[o][d] / shortest_time.max()
                demand[o][d] = ((G.node[origins[o]][weight] * G.node[destinations[d]][weight]) * np.exp(-1 * dist_decay * normalized_dist))

    demand = ((demand / demand.max()) * maxtrips)
    demand = np.ceil(demand).astype(int)

def reflect_roads(G):
    warnings.warn("WARNING! This function is deprecated and will be removed in \
    future releases of GOSTnets. Consider using add_missing_reflected_edges \
    instead", DeprecationWarning)
    #### Function for ensuring bi-directionality of roads ####
    # REQUIRED: G - a graph containing one or more nodes and one or more edges
    # -------------------------------------------------------------------------#
    G_copy = G.copy()

    new_edge_bucket = []

    edge_id = 0

    for u,v, data in G_copy.edges(data = True):
        data['edge_id'] = edge_id
        new_edge = (v,u,data)
        new_edge_bucket.append(new_edge)
        edge_id += 1

    G_copy.add_edges_from(new_edge_bucket)

    return G_copy

def unbundle_geometry(c):
    #### Function for unbundling complex geometric objects ####
    # REQUIRED: any object. This helper function is usually applied in lambda
    #           format against a pandas / geopandas dataframe
    # RETURNS:  an unbundled geometry value that can be plotted.
    # NOTE:     shapely MultiLineString objects quickly get complicated. They
    #           may not show up when you plot them in QGIS. This function aims
    #           to make a .csv 'plottable'
    # -------------------------------------------------------------------------#

    if type(c) == list:
        objs = []
        for i in c:
            if type(i) == str:
                J = loads(i)
                if type(J) == LineString:
                    objs.append(J)
                if type(J) == MultiLineString:
                    for j in J:
                        objs.append(j)
            elif type(i) == MultiLineString:
                for j in i:
                    objs.append(j)
            elif type(i) == LineString:
                objs.append(i)
            else:
                pass
            mls = MultiLineString(objs)
            ls = linemerge(mls)
        return ls
    else:
        return loads(c)

def save(G, savename, wpath):

    ### function used to save a graph object in a variety of handy formats ###
    # REQUIRED:     G - a graph object
    #               savename - the filename, WITHOUT extension
    #               wpath - the write path for where the user wants the files saved
    # -------------------------------------------------------------------------#

    new_node_gdf = node_gdf_from_graph(G)
    new_node_gdf.to_csv(os.path.join(wpath, '%s_nodes.csv' % savename))

    new_edge_gdf = edge_gdf_from_graph(G)
    new_edge_gdf.to_csv(os.path.join(wpath, '%s_edges.csv' % savename))

    nx.write_gpickle(G, os.path.join(wpath, '%s.pickle' % savename))


def add_missing_reflected_edges(G):

    ### function for adding any missing reflected edges - makes all edges
    #   bidirectional. This is essential for routing with simplified graphs ###
    # REQUIRED:     G - a graph object
    # -------------------------------------------------------------------------#

    unique_edges = []
    missing_edges = []

    for u, v, data in G.edges(data = True):
        unique_edges.append((u,v))
    for u, v, data in G.edges(data = True):
        if (v, u) not in unique_edges:
            missing_edges.append((v,u,data))
    G2 = G.copy()
    G2.add_edges_from(missing_edges)
    print(G2.number_of_edges())
    return G2

def remove_duplicate_edges(G):

    ### function for adding any deleting duplicated edges - where there is more
    #   than one edge connecting a node pair. USE WITH CAUTION - will change both
    #   topological relationships and node maps                              ###
    # REQUIRED:     G - a graph object
    # -------------------------------------------------------------------------#

    G2 = G.copy()
    uniques = []
    deletes = []
    for u, v, data in G2.edges(data = True):
        if (u,v) not in uniques:
            uniques.append((u,v))
        else:
            t = G2.number_of_edges(u, v)
            lengths = []
            for i in range(0,t):
                lengths.append(G2.edges[u,v,i]['length'])
            if max(lengths) / min(lengths) >= 1.5:
                pass
            else:
                deletes.append((u,v))
    for d in deletes:
        G2.remove_edge(d[0],d[1])
    print(G2.number_of_edges())
    return G2

def convert_to_MultiDiGraph(G):
    ### takes any graph object, loads it into a MultiDiGraph type Networkx object.
    # REQUIRED:     G - a graph object
    # -------------------------------------------------------------------------#

    a = nx.MultiDiGraph()

    node_bunch = []
    for u, data in G.nodes(data = True):
        node_bunch.append((u,data))

    a.add_nodes_from(node_bunch)

    edge_bunch = []
    for u, v, data in G.edges(data = True):
        data['Wkt'] = str(data['Wkt'])
        edge_bunch.append((u,v,data))

    a.add_edges_from(edge_bunch)
    print(a.number_of_edges())
    return a


#### NETWORK SIMPLIFICATION ####

def simplify_junctions(G, measure_crs, in_crs = {'init': 'epsg:4326'}, thresh = 25):


    ### simplifies topology of networks by simplifying node clusters into single
    # nodes
    # REQUIRED:     G - a graph object
    #               measure_crs - the crs to make the measurements inself.
    # OPTIONAL:    in_crs - the current crs of the graph's geometry properties.
    #               by default, assumes WGS 84 (epsg 4326)
    #               thresh - the threshold distance in which to simplify junctions.
    #               by default, assumes 25 metres
    # -------------------------------------------------------------------------#

    G2 = G.copy()

    gdfnodes = node_gdf_from_graph(G2)
    gdfnodes_proj_buffer = gdfnodes.to_crs(measure_crs)
    gdfnodes_proj_buffer = gdfnodes_proj_buffer.buffer(thresh)
    juncs_gdf = gpd.GeoDataFrame(pd.DataFrame({'geometry':unary_union(gdfnodes_proj_buffer)}), crs = measure_crs, geometry = 'geometry')
    juncs_gdf['area'] = juncs_gdf.area

    juncs_gdf_2 = juncs_gdf.copy()
    juncs_gdf_2 = juncs_gdf_2.loc[juncs_gdf_2.area > int(juncs_gdf.area.min() + 1)]
    juncs_gdf = juncs_gdf_2
    juncs_gdf = juncs_gdf.reset_index()
    juncs_gdf['obj_ID'] = juncs_gdf.index
    juncs_gdf['obj_ID'] = 'new_obj_'+juncs_gdf['obj_ID'].astype(str)

    juncs_gdf_unproj = juncs_gdf.to_crs(in_crs)
    juncs_gdf_unproj['centroid'] = juncs_gdf_unproj.centroid
    juncs_gdf_bound = gpd.sjoin(juncs_gdf_unproj, gdfnodes, how='left', op='intersects', lsuffix='left', rsuffix='right')
    juncs_gdf_bound = juncs_gdf_bound[['obj_ID','centroid','node_ID']]

    node_map = juncs_gdf_bound[['obj_ID','node_ID']]
    node_map = node_map.set_index('node_ID')
    node_dict = node_map['obj_ID'].to_dict()
    nodes_to_be_destroyed = list(node_dict.keys())

    centroid_map = juncs_gdf_bound[['obj_ID','centroid']]
    centroid_map = centroid_map.set_index('obj_ID')
    centroid_dict = centroid_map['centroid'].to_dict()
    new_node_IDs = list(centroid_dict.keys())

    # Add the new centroids of the junction areas as new nodes
    new_nodes = []
    for i in new_node_IDs:
        new_nodes.append((i, {'x':centroid_dict[i].x, 'y':centroid_dict[i].y}))
    G2.add_nodes_from(new_nodes)

    # modify edges - delete those where both u and v are to be removed, edit the others
    edges_to_be_destroyed = []
    new_edges = []

    for u, v, data in G2.edges(data = True):

        if type(data['Wkt']) == LineString:
            l = data['Wkt']
        else:
            l = loads(data['Wkt'])

        line_to_be_edited = l.coords

        if u in nodes_to_be_destroyed and v in nodes_to_be_destroyed:
            if node_dict[u] == node_dict[v]:
                edges_to_be_destroyed.append((u,v))

            else:
                new_ID_u = node_dict[u]
                new_point_u = centroid_dict[new_ID_u]
                new_ID_v = node_dict[v]
                new_point_v = centroid_dict[new_ID_v]

                if len(line_to_be_edited) > 2:
                    data['Wkt'] = LineString([new_point_u, *line_to_be_edited[1:-1], new_point_v])
                else:
                    data['Wkt'] = LineString([new_point_u, new_point_v])
                data['Type'] = 'dual_destruction'

                new_edges.append((new_ID_u,new_ID_v,data))
                edges_to_be_destroyed.append((u,v))

        else:

            if u in nodes_to_be_destroyed:
                new_ID_u = node_dict[u]
                u = new_ID_u

                new_point = centroid_dict[new_ID_u]
                coords = [new_point, *line_to_be_edited[1:]]
                data['Wkt'] = LineString(coords)
                data['Type'] = 'origin_destruction'

                new_edges.append((new_ID_u,v,data))
                edges_to_be_destroyed.append((u,v))

            elif v in nodes_to_be_destroyed:
                new_ID_v = node_dict[v]
                v = new_ID_v

                new_point = centroid_dict[new_ID_v]
                coords = [*line_to_be_edited[:-1], new_point]
                data['Wkt'] = LineString(coords)
                data['Type'] = 'destination_destruction'

                new_edges.append((u,new_ID_v,data))
                edges_to_be_destroyed.append((u,v))

            else:
                data['Type'] = 'legitimate'
                pass

    # remove old edges that connected redundant nodes to each other / edges where geometry needed to be changed
    G2.remove_edges_from(edges_to_be_destroyed)

    # ... and add any corrected / new edges
    G2.add_edges_from(new_edges)

    # remove now redundant nodes
    G2.remove_nodes_from(nodes_to_be_destroyed)

    print(G2.number_of_edges())

    return G2


def custom_simplify(G, strict=True):
    """
    Simplify a graph's topology by removing all nodes that are not intersections
    or dead-ends.

    Create an edge directly between the end points that encapsulate them,
    but retain the geometry of the original edges, saved as attribute in new
    edge.

    Parameters
    ----------
    G : networkx multidigraph
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules
        but have edges with different OSM IDs

    Returns
    -------
    networkx multidigraph
    """

    def get_paths_to_simplify(G, strict=True):

        """
        Create a list of all the paths to be simplified between endpoint nodes.

        The path is ordered from the first endpoint, through the interstitial nodes,
        to the second endpoint. If your street network is in a rural area with many
        interstitial nodes between true edge endpoints, you may want to increase
        your system's recursion limit to avoid recursion errors.

        Parameters
        ----------
        G : networkx multidigraph
        strict : bool
            if False, allow nodes to be end points even if they fail all other rules
            but have edges with different OSM IDs

        Returns
        -------
        paths_to_simplify : list
        """

        # first identify all the nodes that are endpoints
        start_time = time.time()
        endpoints = set([node for node in G.nodes() if is_endpoint(G, node, strict=strict)])

        start_time = time.time()
        paths_to_simplify = []

        # for each endpoint node, look at each of its successor nodes
        for node in endpoints:
            for successor in G.successors(node):
                if successor not in endpoints:
                    # if the successor is not an endpoint, build a path from the
                    # endpoint node to the next endpoint node
                    try:
                        path = build_path(G, successor, endpoints, path=[node, successor])
                        paths_to_simplify.append(path)
                    except RuntimeError:
                        # recursion errors occur if some connected component is a
                        # self-contained ring in which all nodes are not end points.
                        # could also occur in extremely long street segments (eg, in
                        # rural areas) with too many nodes between true endpoints.
                        # handle it by just ignoring that component and letting its
                        # topology remain intact (this should be a rare occurrence)
                        # RuntimeError is what Python <3.5 will throw, Py3.5+ throws
                        # RecursionError but it is a subtype of RuntimeError so it
                        # still gets handled
                        pass

        return paths_to_simplify

    def is_endpoint(G, node, strict=True):
        """
        Return True if the node is a "real" endpoint of an edge in the network, \
        otherwise False. OSM data includes lots of nodes that exist only as points \
        to help streets bend around curves. An end point is a node that either: \
        1) is its own neighbor, ie, it self-loops. \
        2) or, has no incoming edges or no outgoing edges, ie, all its incident \
            edges point inward or all its incident edges point outward. \
        3) or, it does not have exactly two neighbors and degree of 2 or 4. \
        4) or, if strict mode is false, if its edges have different OSM IDs. \

        Parameters
        ----------
        G : networkx multidigraph

        node : int
            the node to examine
        strict : bool
            if False, allow nodes to be end points even if they fail all other rules \
            but have edges with different OSM IDs

        Returns
        -------
        bool

        """
        neighbors = set(list(G.predecessors(node)) + list(G.successors(node)))
        n = len(neighbors)
        d = G.degree(node)

        if node in neighbors:
            # if the node appears in its list of neighbors, it self-loops. this is
            # always an endpoint.
            return 'node in neighbours'

        # if node has no incoming edges or no outgoing edges, it must be an endpoint
        #elif G.out_degree(node)==0 or G.in_degree(node)==0:
            #return 'no in or out'

        elif not (n==2 and (d==2 or d==4)):
            # else, if it does NOT have 2 neighbors AND either 2 or 4 directed
            # edges, it is an endpoint. either it has 1 or 3+ neighbors, in which
            # case it is a dead-end or an intersection of multiple streets or it has
            # 2 neighbors but 3 degree (indicating a change from oneway to twoway)
            # or more than 4 degree (indicating a parallel edge) and thus is an
            # endpoint
            return 'condition 3'

        elif not strict:
            # non-strict mode
            osmids = []

            # add all the edge OSM IDs for incoming edges
            for u in G.predecessors(node):
                for key in G[u][node]:
                    osmids.append(G.edges[u, node, key]['osmid'])

            # add all the edge OSM IDs for outgoing edges
            for v in G.successors(node):
                for key in G[node][v]:
                    osmids.append(G.edges[node, v, key]['osmid'])

            # if there is more than 1 OSM ID in the list of edge OSM IDs then it is
            # an endpoint, if not, it isn't
            return len(set(osmids)) > 1

        else:
            # if none of the preceding rules returned true, then it is not an endpoint
            return False

    def build_path(G, node, endpoints, path):
        """
        Recursively build a path of nodes until you hit an endpoint node.

        Parameters
        ----------
        G : networkx multidigraph
        node : int
            the current node to start from
        endpoints : set
            the set of all nodes in the graph that are endpoints
        path : list
            the list of nodes in order in the path so far

        Returns
        -------
        paths_to_simplify : list
        """
        # for each successor in the passed-in node
        for successor in G.successors(node):
            if successor not in path:
                # if this successor is already in the path, ignore it, otherwise add
                # it to the path
                path.append(successor)
                if successor not in endpoints:
                    # if this successor is not an endpoint, recursively call
                    # build_path until you find an endpoint
                    path = build_path(G, successor, endpoints, path)
                else:
                    # if this successor is an endpoint, we've completed the path,
                    # so return it
                    return path

        if (path[-1] not in endpoints) and (path[0] in G.successors(path[-1])):
            # if the end of the path is not actually an endpoint and the path's
            # first node is a successor of the path's final node, then this is
            # actually a self loop, so add path's first node to end of path to
            # close it
            path.append(path[0])

        return path

    ## MAIN PROCESS FOR CUSTOM SIMPLIFY ##

    G = G.copy()

    if type(G) != nx.classes.multidigraph.MultiDiGraph:
        G = ConvertToMultiDiGraph(G)
    initial_node_count = len(list(G.nodes()))
    initial_edge_count = len(list(G.edges()))
    all_nodes_to_remove = []
    all_edges_to_add = []

    # construct a list of all the paths that need to be simplified
    paths = get_paths_to_simplify(G, strict=strict)

    start_time = time.time()
    for path in paths:

        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        edge_attributes = {}
        for u, v in zip(path[:-1], path[1:]):

            # there shouldn't be multiple edges between interstitial nodes
            if not G.number_of_edges(u, v) == 1:
                pass
            # the only element in this list as long as above check is True
            # (MultiGraphs use keys (the 0 here), indexed with ints from 0 and
            # up)
            edge = G.edges[u, v, 0]
            for key in edge:
                if key in edge_attributes:
                    # if this key already exists in the dict, append it to the
                    # value list
                    edge_attributes[key].append(edge[key])
                else:
                    # if this key doesn't already exist, set the value to a list
                    # containing the one value
                    edge_attributes[key] = [edge[key]]

        for key in edge_attributes:
            # don't touch the length attribute, we'll sum it at the end
            if key == 'Wkt':
                edge_attributes['Wkt'] = list(edge_attributes['Wkt'])
            elif key != 'length' and key != 'Wkt':      # if len(set(edge_attributes[key])) == 1 and not key == 'length':
                # if there's only 1 unique value in this attribute list,
                # consolidate it to the single value (the zero-th)
                edge_attributes[key] = edge_attributes[key][0]
            elif not key == 'length':
                # otherwise, if there are multiple values, keep one of each value
                edge_attributes[key] = list(set(edge_attributes[key]))

        # construct the geometry and sum the lengths of the segments
        edge_attributes['geometry'] = LineString([Point((G.nodes[node]['x'], G.nodes[node]['y'])) for node in path])
        edge_attributes['length'] = sum(edge_attributes['length'])

        # add the nodes and edges to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append({'origin':path[0],
                                 'destination':path[-1],
                                 'attr_dict':edge_attributes})

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge['origin'], edge['destination'], **edge['attr_dict'])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    msg = 'Simplified graph (from {:,} to {:,} nodes and from {:,} to {:,} edges) in {:,.2f} seconds'
    return G

def salt_long_lines(G, source, target, thresh = 5000, factor = 1):

    ### adds in new nodes to edges greater than a given length ###
    # REQUIRED:     G - a graph object
    #               source - crs object in format 'epsg:4326'
    #               target - crs object in format 'epsg:32638'
    # OPTIONAL:    thresh - distance in metres after which to break edges.
    #              factor - edge lengths can be returned in units other than
    #              metres by specifying a numerical multiplication factor
    # -------------------------------------------------------------------------#

    def cut(line, distance):
        # Cuts a line in two at a distance from its starting point
        if distance <= 0.0 or distance >= line.length:
            return [LineString(line)]

        coords = list(line.coords)
        for i, p in enumerate(coords):
            pd = line.project(Point(p))
            if pd == distance:
                return [LineString(coords[:i+1]),LineString(coords[i:])]
            if pd > distance:
                cp = line.interpolate(distance)
                return [LineString(coords[:i] + [(cp.x, cp.y)]),LineString([(cp.x, cp.y)] + coords[i:])]

    G2 = G.copy()

    # define transforms for exchanging between source and target projections
    project_WGS_UTM = partial(
                pyproj.transform,
                pyproj.Proj(init=source),
                pyproj.Proj(init=target))

    project_UTM_WGS = partial(
                pyproj.transform,
                pyproj.Proj(init=target),
                pyproj.Proj(init=source))

    long_edges, long_edge_IDs, unique_long_edges, new_nodes, new_edges = [], [], [], [], []

    # Identify long edges
    for u, v, data in G2.edges(data = True):

        # load geometry
        if type(data['Wkt']) == str:
            WGS_geom = loads(data['Wkt'])
        else:
            WGS_geom = unbundle_geometry(data['Wkt'])
        UTM_geom = transform(project_WGS_UTM, WGS_geom)

        # test geomtry length
        if UTM_geom.length > thresh:
            long_edges.append((u, v, data))
            long_edge_IDs.append((u,v))
            if (v, u) in long_edge_IDs:
                pass
            else:
                unique_long_edges.append((u, v, data))

    print('Identified %d unique edge(s) longer than %d. \nBeginning new node creation...' % (len(unique_long_edges), thresh))

    # iterate through one long edge for each bidirectional long edge pair

    j,o = 1, 0

    for u, v, data in unique_long_edges:

        # load geometry of long edge
        if type(data['Wkt']) == str:
            WGS_geom = loads(data['Wkt'])
        else:
            WGS_geom = unbundle_geometry(data['Wkt'])

        if WGS_geom.type == 'MultiLineString':
            WGS_geom = linemerge(WGS_geom)

        UTM_geom = transform(project_WGS_UTM, WGS_geom)

        # flip u and v if Linestring running from v to u, coordinate-wise
        u_x_cond = round(WGS_geom.coords[0][0], 6) == round(G.nodes()[u]['x'], 6)
        u_y_cond = round(WGS_geom.coords[0][1], 6) == round(G.nodes()[u]['y'], 6)

        v_x_cond = round(WGS_geom.coords[0][0], 6) == round(G.nodes()[v]['x'], 6)
        v_y_cond = round(WGS_geom.coords[0][1], 6) == round(G.nodes()[v]['y'], 6)

        if u_x_cond and u_y_cond:
            pass
        elif v_x_cond and v_y_cond:
            u, v = v, u
        else:
            print('ERROR! FUCKED!')

        # calculate number of new nodes to add along length
        number_of_new_points = UTM_geom.length / thresh

        # for each new node
        for i in range(0, int(number_of_new_points+1)):

            ## GENERATE NEW NODES ##

            cur_dist = (thresh * (i+1))

            # generate new geometry along line
            new_point = UTM_geom.interpolate(cur_dist)

            new_point_WGS = transform(project_UTM_WGS, new_point)

            node_data = {'geometry': new_point_WGS,
                        'x' : new_point_WGS.x,
                        'y': new_point_WGS.y}

            new_node_ID = str(u)+'_'+str(i+j)+'_'+str(o)

            # generate a new node as long as we aren't the final node
            if i < int(number_of_new_points):
                new_nodes.append((new_node_ID, node_data))

            ## GENERATE NEW EDGES ##

            # define geometry to be cutting (iterative)
            if i == 0:
                geom_to_split = UTM_geom

            else:
                geom_to_split = result[1]

            # cut geometry. result[0] is the section cut off, result[1] is remainder
            result = cut(geom_to_split, (thresh))

            edge_data = {'Wkt' : transform(project_UTM_WGS, result[0]),
                        'osm_id' : data['osm_id'],
                        'length' : (factor * int(result[0].length)),
                        'infra_type' : data['infra_type'],
                        }

            if i == 0:
                prev_node_ID = u

            if i == int(number_of_new_points):
                new_node_ID = v

            # append resulting edges to a list of new edges, bidirectional.
            new_edges.append((prev_node_ID,new_node_ID,edge_data))
            new_edges.append((new_node_ID,prev_node_ID,edge_data))

            o += 1

            prev_node_ID = new_node_ID

        j+=1

    # add new nodes and edges
    G2.add_nodes_from(new_nodes)
    G2.add_edges_from(new_edges)

    # remove the too-long edges
    for d in long_edges:
        G2.remove_edge(d[0],d[1])

    print('%d new edges added and %d removed to bring total edges to %d' % (len(new_edges),len(long_edges),G2.number_of_edges()))
    print('%d new nodes added to bring total nodes to %d' % (len(new_nodes),G2.number_of_nodes()))

    return G2

def pandana_snap(G, point_gdf, source_crs = 'epsg:4326', target_crs = 'epsg:4326', add_dist_to_node_col = False):

    import networkx as nx
    import geopandas as gpd
    from shapely.geometry import Point
    from scipy import spatial
    from functools import partial
    import pyproj
    from shapely.ops import transform

    in_df = point_gdf.copy()
    node_gdf = node_gdf_from_graph(G)

    if add_dist_to_node_col is True:

        project_WGS_UTM = partial(
                    pyproj.transform,
                    pyproj.Proj(init=source_crs),
                    pyproj.Proj(init=target_crs))

        in_df = point_gdf.copy()
        in_df['Proj_geometry'] = in_df.apply(lambda x: transform(project_WGS_UTM, x.geometry), axis = 1)
        in_df = in_df.set_geometry('Proj_geometry')
        in_df['x'] = in_df.Proj_geometry.x
        in_df['y'] = in_df.Proj_geometry.y

        node_gdf = node_gdf_from_graph(G)
        node_gdf['Proj_geometry'] = node_gdf.apply(lambda x: transform(project_WGS_UTM, x.geometry), axis = 1)
        node_gdf = node_gdf.set_geometry('Proj_geometry')
        node_gdf['x'] = node_gdf.Proj_geometry.x
        node_gdf['y'] = node_gdf.Proj_geometry.y

        G_tree = spatial.KDTree(node_gdf[['x','y']].as_matrix())

        distances, indices = G_tree.query(in_df[['x','y']].as_matrix())

        in_df['NN'] = list(node_gdf['node_ID'].iloc[indices])
        in_df['NN_dist'] = distances
        in_df = in_df.drop(['x','y'], axis = 1)

    else:
        in_df['x'] = in_df.geometry.x
        in_df['y'] = in_df.geometry.y
        G_tree = spatial.KDTree(node_gdf[['x','y']].as_matrix())

        distances, indices = G_tree.query(in_df[['x','y']].as_matrix())

        in_df['NN'] = list(node_gdf['node_ID'].iloc[indices])

    return in_df
