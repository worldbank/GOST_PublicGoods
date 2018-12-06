import peartree as pt
print('peartree version: %s ' % pt.__version__)
import networkx as nx
print('networkx version: %s ' % nx.__version__)
import matplotlib as mpl
print('matplotlib version: %s ' % mpl.__version__)
import osmnx as ox
print('osmnx version: %s ' % ox.__version__)
import os, sys
import pandas as pd, geopandas as gpd
from shapely.geometry import Point

def node_gdf_from_graph(G, crs = {'init' :'epsg:4326'}, attr_list = None):

    #### Function for generating GeoDataFrame from Graph ####
    # REQUIRED: a graph object G
    # OPTIONAL: crs - projection of format {'init' :'epsg:4326'}. Defaults to WGS84.
    #           note: defining crs of the data, does NOT reproject to this crs
    #           attr_list: list of the keys which you want to be moved over to the GeoDataFrame, if not all
    # RETURNS: a geodataframe of the node objects in the graph
    # -------------------------------------------------------------------------#

    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point

    nodes = []
    keys = []

    if attr_list is None:
        for u, data in G.nodes(data = True):
            keys.append(list(data.keys()))
        flatten = lambda l: [item for sublist in l for item in sublist]
        keys = list(set(flatten(keys)))
        attr_list = keys

    z = 0
    for u, data in G.nodes(data=True):
        if 'geometry' not in attr_list and 'x' in attr_list and 'y' in attr_list :
            new_column_info = {
            'node_ID': u,
            'geometry': Point(data['x'], data['y']),
            'x': data['x'],
            'y': data['y']}
        else:
            new_column_info = {
            'node_ID': u,
            'geometry': Point(u[0],u[1]),
            'x':u[0],
            'y':u[1]}
        for i in attr_list:
            try:
                new_column_info[i] = data[i]
            except:
                pass

        nodes.append(new_column_info)
        z += 1

    nodes_df = pd.DataFrame(nodes)
    nodes_df = nodes_df[['node_ID',*attr_list,'geometry']]
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

    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import LineString

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

    nn = []

    for i, row in points.iterrows():

        pointobj = points[geomcol].loc[i]

        point = (pointobj.x, pointobj.y)

        nearest_nodes = get_nearest_nodes(node_df_G1, point, connection_threshold = connection_threshold)

        nn.append(nearest_nodes.end_node.loc[nearest_nodes.length.idxmin()])

    points['Nearest_node'] = nn

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

    import networkx as nx
    import geopandas as gpd

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

    import networkx as nx
    import geopandas as gpd

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

    from collections import Counter

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

    from shapely.geometry import LineString

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

def convert_network_to_time(G, distance_tag, graph_type = 'drive', speed_dict = None, walk_speed = 4.5):

    #### Function for adding a time value to edge dictionaries ####
    # REQUIRED: G - a graph containing one or more nodes
    #           graph_type - flags network type
    #           distance_tag - the key in the dictionary for the field currently containing a distance in meters
    # OPTIONAL: speed_dict - speed dictionary to use. If not supplied, reverts to defaults
    #           walk_speed - specify a walkspeed in km/h
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

        orig_len = data[distance_tag]

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
            highwayclass = data['highway']

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

    import numpy as np
    import peartree.graph as ptg
    from shapely.geometry import LineString

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
    import numpy as np
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

    import numpy as np

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
