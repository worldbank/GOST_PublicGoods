# -*- coding: utf-8 -*-
###################################################################################################
# Calculate a criticality score for various linear features arranged in a network configuration
# Yan Deng, Charles Fox and Ben Stewart, September 2017
# Purpose: determine the criticality score for each linear feature in the network dataset
###################################################################################################

### Setup
##Import necessary Libraries
import os, sys, inspect, logging
import networkx as nx
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.geometry.base
import shapely.wkt
import copy
import time
import datetime
from datetime import datetime
import warnings

## Remove warnings
warnings.filterwarnings("ignore")

## append to system path the module path from TU Delft. In the GitHub version, this sits one level up in PCS/Criticality
# If designing your own implementation, ensure to add the location of the module to the system path (modify variable 'module_path')
module_path = os.path.join(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe() ))[0])))
module_path = os.path.join(os.path.split(module_path)[0], "PCS/Criticality")
if module_path not in sys.path:
    sys.path.append(module_path)

## Import modules developed by TU Delft team
from network_lib import network_prep as net_p
from network_lib import od_prep as od_p

## Control - verbose prints more to console at runtime, helpful for debugging. Similarly, dump will dump useful files to .csv midway through process
verbose = 1
dump = 1

### Main Function
def main(adminIsPoint = False):

    ## Define filepath
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    path = os.path.split(path)[0]

    ## Define dash. This .xlsm includes settings for the criticality script
    dash = os.path.join(path,r'dashboard.xlsm')
    ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)

    ## Define operative district. Note, this parameter can be anything - it is the sub folder in input, runtime where files are drawn from
    district = ctrl['Weight'].loc['DISTRICT']

    # ensure folders exist
    runtime = os.path.join(path, r'PCS\Criticality\runtime\%s\\' % district)

    ## Add logging
    logging.basicConfig(filename = os.path.join(runtime, "PCS_Criticality.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
    logging.info("Starting Criticality Process")
    print "Running: Criticality Analysis on %s. Do not interrupt" % district

    ## Path Settings
    # outputs
    outpath = os.path.join(path, 'Outputs', '%s' % district)

    for d in [outpath, runtime]:
        if not os.path.isdir(d):
            os.mkdir(d)

    ## Input file setting

    # location of OD
    OD_IN = os.path.join(path, 'PCS\Criticality\Input', '%s' % district)

    # location of administrative boundaries file
    DATA_IN = os.path.join(path, 'PCS\Criticality\Data_Layers')
    inAdmin = os.path.join(DATA_IN,'Poverty_Communes_2009.shp')

    # road network import. Must be a .csv including geometry information of roads.
    inNetworkFile = os.path.join(OD_IN, 'Network.csv')

    # set WGS 84 coordinate reference system
    crs_in = {'init': 'epsg:4326'}

    # ensure folders exist
    for d in [outpath, runtime, OD_IN]:
        if not os.path.isdir(d):
            os.mkdir(d)

    # error checking - Check input data existence
    for curFile in [dash, inNetworkFile, inAdmin,DATA_IN,OD_IN]:
        if not os.path.exists(curFile):
            logging.error("No input found: %s" % curFile)
            raise ValueError("No input found: %s" % curFile)

    # import input dataframes - road network and control dashboard
    inNetwork = pd.read_csv(inNetworkFile)
    ctrldf = pd.read_excel(dash, sheetname = "CRITICALITY", index_col = 'COL_ID')

    #Inputs

    # setting network shapefile location
    network = os.path.join(runtime,'Network.shp')

    ## Network Preparation
    # set default iri value as the mean iri of roads for which iri exists.
    fillvalue = inNetwork['iri_med'].mean()

    # fill iri value where missing
    inNetwork['TC_iri_med'] = inNetwork['iri_med'].fillna(fillvalue)

    # set cost of traversing segment according to length and IRI, per settings in the excel dashboard
    inNetwork['total_cost'] = inNetwork['length'] * (ctrldf['Base_cost_km'][0] + (ctrldf['IRI_Coeff'][0] * inNetwork['TC_iri_med']))

    # convert the pandas DataFrame to a GeoDataFrame
    ginNetwork = gpd.GeoDataFrame(inNetwork,crs = crs_in, geometry = inNetwork['Line_Geometry'].map(shapely.wkt.loads))

    # set up Shapefile of road network

    ginNetwork.to_file(network, driver = 'ESRI Shapefile')
    logging.info("Successfully loaded data")

    # Generate admin boundary centroids
    if not adminIsPoint:
        prepareAdminCentroids(ginNetwork, inAdmin, crs_in, os.path.join(OD_IN, 'adm_centroids.shp'))
        logging.info("Created admin centroids")

    # define function for loading origin files into a dictionary. Paramters controlled from dashboard excel
    def makeOrigin(n, ctrldf):
        origindict = {
            'name': ctrldf['OName'][n],
            'file': os.path.join(path,'PCS','Criticality','Input',district,'%s.shp'% ctrldf['OName'][n]),
            'scalar_column': ctrldf['OScalar'][n]
            }
        return origindict

    # define function for loading destination files into a dictionary. Paramters controlled from dashboard excel
    def makeDestination(n, ctrldf):
        destdict = {
            'name': ctrldf['DName'][n],
            'file': os.path.join(path,'PCS','Criticality','Input',district,'%s.shp'% ctrldf['DName'][n]),
            'penalty': ctrldf['DPenalty'][n],
            'importance':ctrldf['DImportance'][n],
            'annual':ctrldf['DAnnual'][n],
            'scalar_column': ctrldf['DScalar'][n]
            }
        return destdict

    # load origins and destinations into dictionary, create dictionaries of each set
    origin_1, origin_2, origin_3, origin_4, origin_5 = makeOrigin(0, ctrldf), makeOrigin(1, ctrldf), makeOrigin(2, ctrldf), makeOrigin(3, ctrldf), makeOrigin(4, ctrldf)
    originlist = {
        '%s' % ctrldf['OName'][0]: origin_1,
        '%s' % ctrldf['OName'][1]: origin_2,
        '%s' % ctrldf['OName'][2]: origin_3,
        '%s' % ctrldf['OName'][3]: origin_4,
        '%s' % ctrldf['OName'][4]: origin_5,
        }
    destination_1, destination_2, destination_3, destination_4, destination_5 = makeDestination(0, ctrldf), makeDestination(1, ctrldf), makeDestination(2, ctrldf), makeDestination(3, ctrldf), makeDestination(4, ctrldf)
    destinationlist = {
        '%s' % ctrldf['DName'][0]: destination_1,
        '%s' % ctrldf['DName'][1]: destination_2,
        '%s' % ctrldf['DName'][2]: destination_3,
        '%s' % ctrldf['DName'][3]: destination_4,
        '%s' % ctrldf['DName'][4]: destination_5,
        }
    logging.info("Opened origins and destinations")

    # Prepation of network via TU Delft code
    gdf_points, gdf_node_pos, gdf = net_p.prepare_centroids_network(origin_1['file'], network)

    gdf.to_csv(os.path.join(r'C:\Users\charl\Documents\GitHub\Criticality\PCS\Criticality\Runtime\[district_1]','gdf.csv'))
    gdf_node_pos.to_csv(os.path.join(r'C:\Users\charl\Documents\GitHub\Criticality\PCS\Criticality\Runtime\[district_1]','gdf_node_pos.csv'))

        # Create Networkx MultiGraph object from the GeoDataFrame
    G = net_p.gdf_to_simplified_multidigraph(gdf_node_pos, gdf, simplify=False)

    # Change the MultiGraph object to Graph object to reduce computation cost
    G_tograph = net_p.multigraph_to_graph(G)
    logging.info('Loaded road network: number of disconnected components is: %d' % nx.number_connected_components(G_tograph))

    # Observe the properties of the Graph object
    nx.info(G_tograph)

    # Take only the largest subgraph with all connected links
    len_old = 0
    for g in nx.connected_component_subgraphs(G_tograph):
        if len(list(g.edges())) > len_old:
            G1 = g
            len_old = len(list(g.edges()))
    G_sub = G1.copy()

    nx.info(G_sub)

    # Save the simplified transport network into a GeoDataFrame
    gdf_sub = net_p.graph_to_df(G_sub)
    blank, gdf_node_pos2, gdf_new = net_p.prepare_newOD(origin_1['file'], gdf_sub)

    #Road Network Graph prep
    G2_multi = net_p.gdf_to_simplified_multidigraph(gdf_node_pos2, gdf_new, simplify=False)

    # Dump files to runtime if dump = 1
    Filedump(gdf_new, 'Road_Lines', runtime)
    Filedump(gdf_node_pos2,'Road_Nodes', runtime)
    G2 = net_p.multigraph_to_graph(G2_multi)
    gdf2 = net_p.graph_to_df(G2)
    nLink = len(G2.edges())

    # open empty lists
    Outputs, cost_list, iso_list = [], [], []

    ## Run the calculateOD function for each combination of origins and destinations specified in the control excel
    # append all outputs to the Outputs, cost_list and iso_list objects just created
    for z in ctrldf.index:
        if (((ctrldf['ComboO'][z]) != 0) & ((ctrldf['ComboD'][z]) != 0) & (pd.notnull(ctrldf['ComboO'][z])) & (pd.notnull(ctrldf['ComboO'][z]))):
            Q = int(ctrldf['ComboNumber'][z])
            logging.info('Computing | combination %s as origin and %s as destination ' % (ctrldf['ComboO'][z],ctrldf['ComboD'][z]))
            xx = calculateOD(originlist['%s' % ctrldf['ComboO'][z]], destinationlist['%s' % ctrldf['ComboD'][z]], Q, gdf_sub, G2, nLink, gdf2, runtime, ctrldf)
            Outputs.append(xx)
            cost_list.append("Social_Cost_%s" % Q)
            iso_list.append("Isolated_Trips_%s" % Q)

    # drop unneccessary columns
    Output = inNetwork.drop(["geometry",'TC_iri_med','total_cost'],axis = 1)

    # for each object in the Outputs list:
    for o_d_calc in range(0,len(Outputs)):

        # Merge the objects together. This creates multiple columns showing each scenario
        Output = Output.merge(Outputs[o_d_calc]['summary'],how = 'left', on = 'ID')

    # sum across the relevant columns - the 'Social_Cost' columns generated above in calculateOD for each O-D file combo
    Output['Cost_total'] = Output[cost_list].sum(axis = 1)

    # sum across the relevant columns - the 'Isolated_Trips' columns generated above in calculateOD for each O-D file combo
    Output['Iso_total']  = Output[iso_list].sum(axis = 1)

    # Generate an overall criticality score for each road based on user input weights between isolated trips and disrupted trips
    Output['CRIT_SCORE'] = (
        ctrldf['Disrupt_Weight'][0] * Output['Cost_total'] +
        ctrldf['Isolate_Weight'][0] * Output['Iso_total']
    )

    # normalize for each road
    Output['CRIT_SCORE'] = ((Output['CRIT_SCORE'] - Output['CRIT_SCORE'].min()) / (Output['CRIT_SCORE'].max() - Output['CRIT_SCORE'].min()))
    logging.info("Calculated PCS Criticality")
    FileOut(Output,'criticality_output', outpath)

def prepareAdminCentroids(ginNetwork, inAdmin, crs_in, outputFile):
    ## function to clip out from a larger admin network file the bit for the road network in question, and then generate centroids

    # transforms the road network into a single object
    SingleNetworkObj = pd.DataFrame([str(ginNetwork.unary_union)], columns = ['geom'])

    # loads this object in as a GeoDataFrame
    gSingleNetworkObj = gpd.GeoDataFrame(SingleNetworkObj,crs = crs_in, geometry = SingleNetworkObj['geom'].map(shapely.wkt.loads))

    # read administrative boundary file in as a GeoDataFrame
    ginAdmin = gpd.read_file(inAdmin)

    # Ensure projection is common to WGS 84
    ginAdmin = ginAdmin.to_crs(crs_in)

    if 'Pop' not in ginAdmin.columns:
        Vprint('ERROR: Admin boundary file requires a numerical column labelled "Pop" denoting district populations!')
        raise ValueError("ERROR: Admin boundary file requires a numerical column labelled 'Pop' denoting district populations!")

    # drop all polygons which do not intersect the single road network object
    ginAdmin['ID'] = ginAdmin.index
    SelectedAdmins = gpd.sjoin(gSingleNetworkObj, ginAdmin, how="inner",op='intersects')
    SelectedAdmins = SelectedAdmins.drop(['geom','geometry'], axis = 1)
    ginAdmin = ginAdmin.loc[ginAdmin['ID'].isin(SelectedAdmins['ID']) == True]

    # transform geoemetry column to centroid
    ginAdmin['geometry'] = ginAdmin['geometry'].centroid

    # print file to shapefile with name outputFile
    ginAdmin.to_file(outputFile, driver = 'ESRI Shapefile')

#Utility Functions
def Filedump(df, name, runtime):
    # Dumps file to runtime folder if dump = 1
    if dump == 1:
        df.to_csv(os.path.join(runtime, r'%s.csv' % name), encoding = 'utf-8')

def FileOut(df, name, outpath):
    # Dumps file to outpath folder
    df.to_csv(os.path.join(outpath, r'%s.csv' % name), index = False, encoding = 'utf-8')
    try:
        df = df.drop('Line_Geometry', axis = 1)
        df.to_excel(os.path.join(outpath, r'%s.xlsx' % name), index = False, encoding = 'utf-8')
    except:
        pass

def Vprint(s):
    # prints statement with timestamp if verbose = 1
    if verbose == 1:
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print '\n%s -- %s' % (st, s)

##########Criticality Functions #########

def PrepSet(point_set, gdf_sub):
    '''
    Prepares a small df of a given origin / destination set, expressed as 'item : Nearest Node ID'
    '''
    Prepared_point_set, gdf_node_pos2, gdf_new = net_p.prepare_newOD(point_set['file'], gdf_sub)
    Prepared_point_set = Prepared_point_set['Node']
    return Prepared_point_set

def Pathfinder(origin, destination, G2):
    '''
    Find the shortest path for each OD to minimize the total travel cost.
    output:
    1) baseCost ($): total travel cost between all OD pairs;
    2) basePath: the shortest path between all OD pairs;
    3) baseLength: total travel distance between all OD pairs.
    '''
    basePath = [[ [] for d in range(len(destination))] for d in range(len(origin))]
    baseCost = np.zeros((len(origin),len(destination)))
    baseLength = np.zeros((len(origin),len(destination)))
    for o in range(len(origin)):
        for d in range(len(destination)):
            basePath[o][d] = nx.dijkstra_path(G2,origin[o],destination[d], weight = 'total_cost')
            baseCost[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'total_cost')
            baseLength[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'length')
    return basePath, baseCost, baseLength

def BreakEdge(origin, destination, penalty, baseCost, name, G2, runtime, nLink, gdf2):
    '''
    Fundamental to the criticaility analysis: remove each link by changing its traverse cost to an extremely large value,
    calculate the shortest path, then re-set to the original cost of traversing the link.
    output:
    1.) diff: a 3D matrix describing the change in cost for each 2D O-D matrix when a given link is broken
    2.) iso: a 3D matrix describing the cost of isolated journeys for each 2D O-D matrix when a given link is broken
    '''
    to_allNode = []
    G = copy.deepcopy(G2)

    # generate an empty 3D matrix - number of road links x origins x destinations
    cost_disrupt = np.zeros((nLink, len(origin),len(destination)))
    ts = time.time()

    # iterate through each road in the network GeoDataFrame:
    for road in range(len(gdf2)):
        ts1 = time.time()
        # w is the original traverse cost, stored as variable for later
        w = G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost']

        # link traverse cost raised to 1e10
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = 1e10
        for o in range(len(origin)):
            for d in range(len(destination)):
                cost_disrupt[road][o][d] = nx.dijkstra_path_length(G,origin[o],destination[d], weight = 'total_cost')
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = w # Resetting traverse cost to original (w)
        logging.info("Computation completed for road: %s of %s. Compute time: %s seconds" % (road+1, len(gdf2)+1, (time.time() - ts1)))

    # disruption cost matrix for a given road, expressed in money
    diff = (cost_disrupt - baseCost)

    # matrix recording isolated journeys by number (each isolated journey = 1)
    iso = np.zeros((nLink, len(origin),len(destination)))

    # change cost of the isolated trips in the OD matrix to the penalty value
    for index,item in np.ndenumerate(cost_disrupt):
        if item >= 1e9:
            diff[index] = penalty
            iso[index] = 1
    return diff, iso

def GenerateDemandFunction(origin,destination,baseLength):
    '''
    Criticality requires the modelling of trips between origins and destinations.
    This function generates demand between each origin and destination according to a gravity model
    The weighting of the origins and destinations is defined in the excel dashboard.
    1.) Output: a demand function of form:
    Trips[o,d] = (Pop[o] * Annual_Number_of_trips * x)
        where sum(x) is 1 and  x = f(norm(Scalar[d]), (exp -(Distance[o,d] * importance))
    '''
    demand = np.zeros((len(origin['P']),len(destination['P'])))
    O_DF = gpd.read_file(origin['file'])
    D_DF = gpd.read_file(destination['file'])
    for o in O_DF.index:

        ## set up destination weighting series for each origin to all destinations
        # dest scalar: how important is each destination? Set in each destination file. Refers to scalar column set in dash excel
        destscalar = pd.Series(D_DF[destination['scalar_column']] / sum(D_DF[destination['scalar_column']]))

        # dist scalar: how far away from o is each destination? Exponential decay model used.
        distscalar = pd.Series(np.exp(-1*(baseLength[o] * (destination['importance']/10))))

        # error catching
        distscalar[distscalar >= 1] = 0

        # destination weighting is a mutliplier of distance and destination importance
        weight = pd.Series(destscalar * distscalar)

        # normalize
        normweight = pd.Series(weight / sum(weight))

        # X is now the dataframe of destination weights from origin o. This changes for each origin.
        X = pd.DataFrame({'x': normweight})

        # for each destination:
        for d in D_DF.index:
            if (o == d) & (origin['name'] == destination['name']):
                pass
            else:
                # generate the value in the trip demand matrix
                demand[o][d] = ((O_DF[origin['scalar_column']][o] * destination['annual']) * X['x'][d])

    # return matrix of all modelled trips - for that particular Origin:Destination combo.
    return demand

def summarise(diff, iso, demand, origin, destination, Q, nLink, gdf2, ctrldf):
    '''
    For each link disrupted, we find the number of trips that cannot be completed due to link disruption: isolate_sumTrip
    Total social cost is defined as the additional cost incurred as jounryes are completed due to the disruption of this link: disrupt_sumCost
    Total number of trips being disrupted: disrupt_sumTrip
    '''
    disrupt_sumCost, disrupt_sumTrip, isolate_sumTrip = np.zeros(nLink),np.zeros(nLink),np.zeros(nLink)

    # works out total costs and total isolated trips by multiplying the by-link information by the demand for each trip
    for i in range((nLink)):
        disrupt_sumCost[i] = np.nansum(np.multiply(diff[i,:,:],demand))
        isolate_sumTrip[i] = np.nansum(np.multiply(iso[i,:,:],demand))

    # make criticality dataframe, 'out', in format 'ID', 'social cost of trips disrupted', 'number of isolated trips'
    criticality = np.column_stack((gdf2['ID'],disrupt_sumCost,isolate_sumTrip))
    criticality[:,0] = criticality[:,0].astype(int)
    out = pd.DataFrame({'ID':criticality[:,0],'Social_Cost_%s' % Q:criticality[:,1],'Isolated_Trips_%s' % Q:criticality[:,2]})

    # Q is the O-D combo number, passed to the function. Prevents overwriting of data as we iterate through O-D combos
    a = out['Social_Cost_%s' % Q]
    b = out['Isolated_Trips_%s' % Q]

    # normalize scores
    out['Social_Cost_score_%s' % Q] = ((a - a.min()) / (a.max()- a.min()))
    out['Isolated_score_%s' % Q] = ((b - b.min()) / (b.max()- b.min()))

    # generate criticality score based on weightings in excel dashboard
    out['CRIT_SCORE_%s' % Q] = (
        ctrldf['Disrupt_Weight'][0] * out['Social_Cost_score_%s' % Q] +
        ctrldf['Isolate_Weight'][0] * out['Isolated_score_%s' % Q]
        )
    out = out[['ID','Social_Cost_%s' % Q,'Isolated_Trips_%s' % Q,'CRIT_SCORE_%s' % Q]]
    return out

def calculateOD(origin, destination, Q, gdf_sub, G2, nLink, gdf2, runtime, ctrldf):
    '''
    the main function - calls all of the lower order function for a given origin:destination pair.
    returns a data dictionary which is assembled in the final steps of the main function process.
    '''
    # binds origins and destinations onto the graph network - snapping points to nearest node
    origin['P'], destination['P'] = PrepSet(origin, gdf_sub), PrepSet(destination, gdf_sub)

    # generates least cost paths, costs and distances across the road network for each origin to each destination
    basePath, baseCost, baseLength = Pathfinder(origin['P'], destination['P'], G2)

    # performs the criticality analysis - iteratively breaking each link, calculating the cost of disrupted and number of isolated trips
    diff, iso = BreakEdge(origin['P'], destination['P'], destination['penalty'], baseCost, destination['name'], G2, runtime, nLink, gdf2)

    # generates a demand function from the data in the origin and destination files.
    demand = GenerateDemandFunction(origin, destination, baseLength)

    # generates summary results
    summary = summarise(diff, iso, demand, origin, destination, Q, nLink, gdf2, ctrldf)

    # returns results dictionary for a given origin:destination file combination
    return({
        'origin': origin['name'],
        'destination': destination['name'],
        'basePath': basePath,
        'baseCost': baseCost,
        'baseLength': baseLength,
        'diff': diff,
        'iso': iso,
        'summary': summary})

main()
