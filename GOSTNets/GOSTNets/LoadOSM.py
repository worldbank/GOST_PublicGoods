####################################################################################################
# Load OSM data into network graph
# Benjamin Stewart and Charles Fox
# Purpose: take an input dataset as a OSM file and return a network object
####################################################################################################

import os, sys, time

import shapely.ops

import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from osgeo import ogr
from rtree import index
from shapely import speedups
from shapely.geometry import LineString, MultiLineString, MultiPoint, Point
from geopy.distance import vincenty
from boltons.iterutils import pairwise
from shapely.wkt import loads,dumps

class OSM_to_network(object):

    def __init__(self, osmFile):
        ''' Generate a networkX object from a osm file
        '''
        self.osmFile = osmFile
        self.roads_raw = self.fetch_roads(osmFile)

    def generateRoadsGDF(self, in_df = None, outFile='', verbose = False):
        if type(in_df) != gpd.geodataframe.GeoDataFrame:
            in_df = self.roads_raw
        roads = self.get_all_intersections(in_df, verboseness = verbose)
        roads['key'] = ['edge_'+str(x+1) for x in range(len(roads))]
        np.arange(1,len(roads)+1,1)

        def get_nodes(x):
            return list(x.geometry.coords)[0],list(x.geometry.coords)[-1]

        nodes = gpd.GeoDataFrame(roads.apply(lambda x: get_nodes(x),axis=1).apply(pd.Series))
        nodes.columns = ['u','v']

        roads['length'] = roads.geometry.apply(lambda x : self.line_length(x))

        #G = ox.gdfs_to_graph(all_nodes,roads)
        roads.rename(columns={'geometry':'Wkt'}, inplace=True)

        roads = pd.concat([roads,nodes],axis=1)

        if outFile != '':
            roads.to_csv(outFile)

        self.roadsGPD = roads

    def filterRoads(self, acceptedRoads = ['primary','primary_link','secondary','secondary_link','motorway','motorway_link','trunk','trunk_link']):
        self.roads_raw = self.roads_raw.loc[self.roads_raw.infra_type.isin(acceptedRoads)]

    def fetch_roads(self, data_path):

        if data_path.split('.')[-1] == 'pbf':

            driver = ogr.GetDriverByName('OSM')
            data = driver.Open(data_path)

            sql_lyr = data.ExecuteSQL("SELECT osm_id,highway FROM lines WHERE highway IS NOT NULL")

            roads=[]
            for feature in sql_lyr:
                if feature.GetField('highway') is not None:
                    osm_id = feature.GetField('osm_id')
                    shapely_geo = loads(feature.geometry().ExportToWkt())
                    if shapely_geo is None:
                        continue
                    highway=feature.GetField('highway')
                    roads.append([osm_id,highway,shapely_geo])

            if len(roads) > 0:
                road_gdf = gpd.GeoDataFrame(roads,columns=['osm_id','infra_type','geometry'],crs={'init': 'epsg:4326'})
                return road_gdf

        elif data_path.split('.')[-1] == 'shp':
            road_gdf = gpd.read_file(data_path)
            return road_gdf

        else:
            print('No roads found')

    def line_length(self, line, ellipsoid='WGS-84'):
        """Length of a line in meters, given in geographic coordinates

        Adapted from https://gis.stackexchange.com/questions/4022/looking-for-a-pythonic-way-to-calculate-the-length-of-a-wkt-linestring#answer-115285

        Arguments:
            line {Shapely LineString} -- a shapely LineString object with WGS-84 coordinates
            ellipsoid {String} -- string name of an ellipsoid that `geopy` understands (see
                http://geopy.readthedocs.io/en/latest/#module-geopy.distance)

        Returns:
            Length of line in meters
        """
        if line.geometryType() == 'MultiLineString':
            return sum(line_length(segment) for segment in line)

        return sum(
                    vincenty(tuple(reversed(a)), tuple(reversed(b)), ellipsoid=ellipsoid).kilometers
                    for a, b in pairwise(line.coords)
        )

    def get_all_intersections(self, shape_input, idx_osm = None, verboseness = False):
        # Initialize Rtree
        idx_inters = index.Index()
        # Load data
        #all_data = dict(zip(list(shape_input.osm_id),list(shape_input.geometry),list(shape_input.infra_type)))
        ### TODO - it shouldn't be necessary to reference the geometry column specifically
        #   ... but here we are
        if idx_osm is None:
            idx_osm = shape_input['geometry'].sindex

        # Find all the intersecting lines to prepare for cutting
        count = 0
        tLength = shape_input.shape[0]
        inters_done = {}
        new_lines = []
        allCounts = []

        for idx, row in shape_input.iterrows():
            key1 = row.osm_id
            line = row.geometry
            infra_type = row.infra_type
            ### TIMING
            if count % 1000 == 0 and verboseness == True:
                print("Processing %s of %s" % (count, tLength))
            count += 1
            #infra_line = line['infra_type']#shape_input.at[shape_input.index[shape_input['osm_id']==key1].tolist()[0],'infra_type']
            ### TIMING
            intersections = shape_input.iloc[list(idx_osm.intersection(line.bounds))]
            intersections = dict(zip(list(intersections.osm_id),list(intersections.geometry)))
            ### TIMING
            # Remove line1
            if key1 in intersections: intersections.pop(key1)
            # Find intersecting lines
            ### TIMING
            for key2,line2 in intersections.items():
                # Check that this intersection has not been recorded already
                if (key1, key2) in inters_done or (key2, key1) in inters_done:
                    continue
                # Record that this intersection was saved
                inters_done[(key1, key2)] = True
                # Get intersection
                if line.intersects(line2):
                    # Get intersection
                    inter = line.intersection(line2)
                    # Save intersecting point
                    if "Point" == inter.type:
                        idx_inters.insert(0, inter.bounds, inter)
                    elif "MultiPoint" == inter.type:
                        for pt in inter:
                            idx_inters.insert(0, pt.bounds, pt)

            # cut lines where necessary and save all new linestrings to a list
            hits = [n.object for n in idx_inters.intersection(line.bounds, objects=True)]

            if len(hits) != 0:
                out = shapely.ops.split(line, MultiPoint(hits))
                new_lines.append([{'geometry': LineString(x), 'osm_id':key1,'infra_type':infra_type} for x in out.geoms])
            else:
                new_lines.append([{'geometry': line, 'osm_id':key1,
                        'infra_type':infra_type}])

        # Create one big list and treat all the cutted lines as unique lines
        flat_list = []
        all_data = {}

        #item for sublist in new_lines for item in sublist
        i = 1
        for sublist in new_lines:
            if sublist is not None:
                for item in sublist:
                    item['id'] = i
                    flat_list.append(item)
                    i += 1
                    all_data[i] = item

        # Transform into geodataframe and add coordinate system
        full_gpd = gpd.GeoDataFrame(flat_list,geometry ='geometry')
        full_gpd.crs = {'init' :'epsg:4326'}
        return(full_gpd)

    def initialReadIn(self, fpath=None):
        if isinstance(fpath, str):
            edges_1 = pd.read_csv(fpath)
            edges_1 = edges_1['Wkt'].apply(lambda x: loads(x))
        elif isinstance(fpath, gpd.GeoDataFrame):
            edges_1 = fpath
        else:
            edges_1 = self.roadsGPD
        edges = edges_1.copy()
        node_bunch = list(set(list(edges['u']) + list(edges['v'])))
        def convert(x):
            u = x.u
            v = x.v
            data = {'Wkt':x.Wkt,
                   'id':x.id,
                   'infra_type':x.infra_type,
                   'osm_id':x.osm_id,
                   'key': x.key,
                   'length':x.length}
            return (u, v, data)

        edge_bunch = edges.apply(lambda x: convert(x), axis = 1).tolist()
        G = nx.MultiDiGraph()
        G.add_nodes_from(node_bunch)
        G.add_edges_from(edge_bunch)
        for u, data in G.nodes(data = True):
            if type(u) == str:
                q = tuple(float(x) for x in u[1:-1].split(','))
            if type(u) == tuple:
                q = u
            data['x'] = q[0]
            data['y'] = q[1]
        G = nx.convert_node_labels_to_integers(G)
        self.network = G
        return G
