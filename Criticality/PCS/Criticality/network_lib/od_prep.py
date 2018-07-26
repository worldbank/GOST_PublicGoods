from __future__ import division
import os
import pandas as pd
import numpy as np
import geopandas as gp

from math import radians, cos, sin, asin, sqrt
from shapely.geometry import LineString, Point

__all__ = ['prepare_centroids_list',
           'calc_distance_centroid',
           'gen_prod_attr',
           'CalcDoublyConstrained',
           'dist_deterrence',
           'district_stats_to_OD_df',
           'all_ods_creation',
           'all_ods_creation_ema',
           'od_ton_to_truck',
           'factors_dict_creation',
           'od_aggregation',
           'od_preparation']

#some codes are obtained from https://github.com/joshchea/python-tdm/blob/master/scripts/CalcDistribution.py
#http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points

def prepare_centroids_list(G2_new_tograph):
    '''
    Extracting list of all centroid nodes' ID from the transport network object

    Parameters
    ------------
    G2_new_tograph: Graph
        Graph Networkx object

    Returns
    ------------
    centroid_nodes: list
        List of all centroid nodes' ID
    '''

    #copy the input object to avoid working directly on the input object
    G = G2_new_tograph.copy()

    #extract subgraph of centroids from the input nodes
    SG=G.subgraph( [n[0] for n in G.node.items() if n[1]['IsCentroid'] == 1 ] )
    SG.nodes(data=True)

    #convert them into list
    centroid_nodes = list(SG.nodes())

    return centroid_nodes

#extract the longitude and latitude from geometry of the shapefile
def _getXY(pt):
    return (pt.x, pt.y)

def _haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    return km

def calc_distance_centroid(centroid_gdf):
    '''
    Calculate euclidean distance between all centroids

    Parameters
    ------------
    centroid_gdf: GeoDataFrame
        GeoDataFrame of centroid which contains Point 'geometry' information

    Returns
    ------------
    distance: DataFrame
        n x n DataFrame of euclidean distance between all centroids
    '''

    #calculate longitude and latitude of points based on shapefile's geometry attribute
    lon,lat = [list(t) for t in zip(*map(_getXY, centroid_gdf['geometry']))]

    #make an arbitrary dataframe to store distance information
    distance = pd.DataFrame({'initiate' : []})

    #calculate the distance between each OD pair
    for i in range(len(lon)):
        d = []
        for j in range(len(lat)):
            d.append(_haversine(lon[i], lat[i], lon[j], lat[j]))
        distance[i] = d
    distance.drop(distance.columns[0], axis=1, inplace=True)

    return distance

def gen_prod_attr(district_stats, prod_driver, attr_driver='Population_x'):
    '''
    Generating single product production and attraction of each district

    Parameters
    ------------
    district_stats: (Geo)DataFrame
        District (Geo)DataFrame which has information about attraction and production
        driver in its columns
    prod_driver: str
        string of factor that will be used as the production driver from the district_stats column
    attr_driver: str
        string of factor that will be used as the attraction driver from the district_stats column

    Returns
    ------------
    production: Series
        Series with centroid id as index and production (in absolute number) as value
    attraction: Series
        Series with centroid id as index and attraction (in relative number) as value
    '''

    #GENERATING TRIP PRODUCTION
    district_stats['trips_production'] = district_stats[prod_driver]
    production = district_stats['trips_production']

    #GENERATING TRIP ATTRACTION
    #if there is no information about the corresponding attraction driver for a particular centroid
    #use mean value of the attraction driver over all centroids
    district_stats[attr_driver] = district_stats[attr_driver].fillna(district_stats[attr_driver].mean())
    #calculate relative attraction of each district
    relative_attr = district_stats[attr_driver] / district_stats[attr_driver].sum()

    #then distribute the production over the relative attraction
    attraction = relative_attr*production.sum()

    return production, attraction

def CalcDoublyConstrained(ProdA, AttrA, F, maxIter = 10):
    '''
    Calculates doubly constrained trip distribution for a given friction factor matrix.

    Parameters
    ------------
    ProdA: Series
        Series of production array from all centroids
    AttrA: Series
        Series of attraction array from all centroids
    F: Series
        Series of friction (penalty) factor for distributing goods from one centroid
        to all other centroids. The factor works in inverse way (e.g. the farther two centroids are located,
        the smaller the friction factor should be).
    maxIter: int
        Maximum iteration of the doubly constrained algorithm

    Returns
    ------------
    Trips1: DataFrame
        DataFrame of OD matrix between all centroids
    '''

    Trips1 = np.zeros((len(ProdA),len(ProdA)))
#     print('Checking production, attraction balancing:')
    sumP = sum(ProdA)
    sumA = sum(AttrA)
    if sumP != sumA:
        AttrA = AttrA*(sumP/sumA)
        AttrT = AttrA.copy()
        ProdT = ProdA.copy()
    else:
        AttrT = AttrA.copy()
        ProdT = ProdA.copy()

    for balIter in range(0, maxIter):
        for i in list(range(0,len(ProdA))):
            Trips1[i,:] = ProdA[i]*AttrA*F[i]/max(0.000001, sum(AttrA * F[i]))

        #Run 2D balancing --->
        ComputedAttractions = Trips1.sum(0)
        ComputedAttractions[ComputedAttractions==0]=1
        AttrA = AttrA*(AttrT/ComputedAttractions)

        ComputedProductions = Trips1.sum(1)
        ComputedProductions[ComputedProductions==0]=1
        ProdA = ProdA*(ProdT/ComputedProductions)

    for i in list(range(0,len(ProdA))):
        c = ProdA[i]*AttrA*F[i]/max(0.000001, sum(AttrA * F[i]))
        Trips1[i,:] = c

    dfc = pd.DataFrame(Trips1)
    Trips1 = dfc.values.tolist()
    return Trips1

def dist_deterrence(distance):
    #Default deterrence function if not specified
    distance = 10000/distance
    for i in list(distance.columns):
        for j in list(distance.index.values):
            if distance[i][j] > 9999999:
                distance[i][j] = 0
    return distance

def district_stats_to_OD_df(gdf_points, prod_driver, attr_driver='Population_x', deterrence=dist_deterrence):
    '''
    Calculating OD flow for one production-attraction pair. Gravity model with euclidean distance
    between centroids are used to calculate the OD flow. Absolute number of the production driver
    is preserved (e.g. goods flow in tonnes).

    Parameters
    ------------
    gdf_points: GeoDataFrame
        geodataframe (Points) of centroids shapefile
    prod_driver: str
        string of factor that will be used as the production driver from the gdf_points column
    attr_driver: str
        string of factor that will be used as the attraction driver from the gdf_points column
    deterrence: function
        deterrence function used to calculate cost value between OD pair

    Returns
    ------------
    OD_matrix: DataFrame
        DataFrame of n x n OD matrix between all centroids
    '''

    #calculate absolute euclidean distance between centroids into DataFrame
    distance = calc_distance_centroid(gdf_points)

    #apply deterrence function to the distance DataFrame
    friction_matrix = deterrence(distance)

    #calcualte production and attraction based on the production driver
    production, attraction = gen_prod_attr(gdf_points, prod_driver, attr_driver)

    #calculate OD_Matrix
    Trips1 = CalcDoublyConstrained(production, attraction, friction_matrix)

    #rename the index and column into nodelist (based on the gdf_points)
    nodelist = list(gdf_points['Node'])
    OD_matrix = pd.DataFrame(Trips1, index=nodelist, columns=nodelist)

    #avoid extremely small number in the flow
    for i, row in OD_matrix.iterrows():
        for j,val in row.iteritems():
            if OD_matrix[i][j] < 0.1:
                OD_matrix[i][j] = 0

    return OD_matrix

def all_ods_creation(gdf_points, prod_lists, attr_driver, dist_deterrence=dist_deterrence):
    '''
    Developing OD matrix based on the socioeconomic data which is available in the centroids GeoDataFrame.
    Gravity function with euclidean distance between centroids is used to calculate the OD flow.
    Absolute number of the production driver is preserved (e.g. goods flow in tonnes).

    Parameters
    ------------
    gdf_points: GeoDataFrame
        geodataframe (Points) of centroids shapefile
    prod_lists: list
        list of string of all production driver. The string comes from the column name of the
        centroids GeoDataFrame
    attr_driver: str
        string of factor that will be used as the attraction driver from the gdf_points column
    dist_deterrence: function
        deterrence function used to calculate penalty between OD pair

    Returns
    ------------
    od_dict: dict
        Dictionary with production drivers as the keys and n x n OD matrix DataFrame between
        all centroids for each production driver as values
    '''

    #create empty dict
    od_dict={}

    #iterate for each production driver
    for prod in prod_lists:
        #calculate the OD matrix
        od = district_stats_to_OD_df(gdf_points,  prod_driver=prod, attr_driver=attr_driver,
                                     deterrence=dist_deterrence)
        #append it to the dict
        od_dict["od_{0}".format(prod)]=od

    return od_dict

def all_ods_creation_ema(gdf_points, prod_lists,attr_driver):

    od_dict={}
    for prod in prod_lists:
        od = district_stats_to_OD_df(gdf_points,  prod_driver=prod, attr_driver=attr_driver)
        od_dict["od_{0}".format(prod)]=(od,prod)

    return od_dict

def _od_ton_to_truck(od, capacity):
    return od / capacity

def _merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def factors_dict_creation(prod_lists):
    #create arbitrary dictionary
    factors_dict={}

    #create scaling factors (for EMA run later)
    factors_scale= [1] * len(prod_lists)

    #enumerate all items in production lists
    for i,prod in enumerate(prod_lists):
        #create new item in dictionary with factor_00, factor_01, etc as keys
        #and production name (e.g. Textile_exp_ton) as values
        factors_dict[prod]=factors_scale[i]

    return factors_dict

def od_aggregation(OD_all_dict, **factors_dict):
    #create empty dictionary
    OD_final_dict={}

    #iterate over all items in original OD
    for key1,val1 in OD_all_dict.iteritems():

        #matching the production value of the OD dict and the factors_dict
        for key2,val2 in factors_dict.iteritems():
            #if it is a match
            if val1[1] == key2:
                #scale the OD flows of that particular product
                OD_final_dict["od_{0}".format(val1[1])]=val1[0]*val2

    #creation of final OD dataframe
    OD_final_df = OD_final_dict[OD_final_dict.keys()[0]]
    for i in range(len(OD_final_dict)-1):
        OD_final_df = OD_final_df +  OD_final_dict[OD_final_dict.keys()[i+1]]

    return OD_final_df

def od_preparation(prod_lists, OD_all_dict, **factors_dict):
    OD_final_df = od_aggregation(OD_all_dict, **factors_dict)
    return OD_final_df