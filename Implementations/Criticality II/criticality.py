
# coding: utf-8
# import libraries
import sys,os
import networkx as nx
import geopandas as gpd
import math
import numpy as np
import time
import pandas as pd
from multiprocessing import Pool
from rasterstats import zonal_stats
import shapely
from shapely.geometry import Point
import subprocess
sys.path.append(os.path.join( '..'))
# append current workin directory to PATH environment variable
sys.path.append(os.getcwd())

# this is a special, local python library called 'utils'. All it does is load file paths from the config.json file
from utils import load_config
import multiprocessing

# import GOSTnets from your location of choice
sys.path.append(r'/home/wb493355/data/criticality II/scripts/')
import GOSTnet as gn

# remove extraneous warning messages with warnings library
import warnings
warnings.filterwarnings("ignore")


def unzip_worldpop(country,base_path,temp_path):

    """Function to unzip the worldpop files for the following islands (not included in the main worldpop file):

        - Fiji
        - Kiribati
        - Marshall Islands
        - Micronesia (Federated States of)
        - Palau
        - Papua New Guinea
        - Samoa
        - Solomon Islands
        - Tonga
        - Vanuatu

    This should work out of the box if *7zip* is installed.

    If not: the easiest way to get this to run, is to move the *7z.exe* into the **scripts** directory. The other option would be to add the *7zip* directory to your environmental variables.

    Args:
        *country* : The ISO-code of the country.

        *base_path* : The base path to the location of all files and scripts.

        *temp_path* : The path to which we temporarily write the unzipped file.

    Returns:
        An unzipped worldpop GeoTIFF file for the specified country.

    """

    archiveman = r'7z' # 7z.exe comes with 7-zip

    """Load file paths"""
    poppath = os.path.join(base_path,'Worldpop')

    islands = ['FJI','KIR','MHL','FSM','PLW','PNG','WSM','SLB','TON','VUT']
    islands_full = ['Fiji','Kiribati','Marshall Islands','Micronesia (Federated States of)','Palau','Papua New Guinea','Samoa','Solomon Islands','Tonga','Vanuatu']

    map_isl_names = dict(zip(islands,islands_full))

    island_full_name = map_isl_names[country]
    archive_name = os.path.join(poppath,'%s 100m Population.7z' % island_full_name)
    if country == 'KIR':
        subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15adj_lzw.tif'])
    elif country == 'TON':
        subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15.tif'])
    else:
        subprocess.check_output([archiveman, 'e','-aos', archive_name, '-o{}'.format(temp_path), 'popmap15adj.tif'])

def RandomDestruction(G, edge_frac, fail_value, seed):
    ### function randomly sets the traverse time on chosen edges to the fail_value, simulating road unavailability
    # edge_frac - proportion of total edges to be disrupted
    # fail_value - large int, traverse time if edge disrupted.
    # seed - random seed. Used to select random values in an ordered way, if necessary

    # make a list of the unique edges by ID.
    # This is NOT the same as a list of all edges - each 'lane' of a road has the same ID but is a separate edge.
    # by selecting based on ID, we prevent one carriageway of a road being disrupted, but the other being open.
    edgeid = []
    for u,v, data in G.edges(data = True):
        edgeid.append(data['edge_id'])
    edgeid = list(set(edgeid))

    # Calculate the number of edges to be destroyed. Always rounds down, if a fraction.
    num_to_destroy = math.floor(len(edgeid) * (edge_frac / 100))

    # set the random seed
    np.random.seed(seed = seed)

    # shuffle the edges list randomly
    np.random.shuffle(edgeid)

    # choose the part of the list up to the number to destroy
    destroy_list = edgeid[:num_to_destroy]

    # create a copy of the graph to adjust it
    G_adj = G.copy()

    # iterate through edge, and set travel time to fail value if in the list of edges to be destroyed
    for u, v, data in G_adj.edges(data = True):
        if data['edge_id'] in destroy_list:
            data['time'] = fail_value

    # function returns the adjusted graph, and the number of edges destroyed
    return G_adj, destroy_list

def Calculate_OD(G, origins, destinations, fail_value, weight = 'time'):
    ### returns numpy O-D array for shortest travel time between each origins and destination.
    ### travel time between o and d is 'answers[o][d]'
    ### incompletable trips return the fail_value
    # G - a graph Object
    # origins - a list of nodes, by ID, that exist in G
    # destinations - a list of nodes, by ID, that exist in G
    # fail_value - must be a large int. Returned when journet between a node pair can't be completed
    # weight - we find shortest paths based on an edge weight - a property of the edges of the graph

    # make an empty numpy array to populate with the travel times
    OD = np.zeros((len(origins), len(destinations)))

    # we iterate through the origins ONLY.
    for o in range(0, len(origins)):

        # the current origin is selected using numerical indexing from the list of origin nodes
        origin = origins[o]

        # then, we calculate the shortest path from this origin to ALL potential destination nodes
        results_dict = nx.single_source_dijkstra_path_length(G, origin, cutoff = None, weight = weight)

        # from this results dictionary, we select just the nodes we need - the destination nodes.
        for d in range(0, len(destinations)):

            # select a single destination from the destinations list using numerical indexing
            destination = destinations[d]

            # populate the OD matrix with the correct result at [o][d]
            if destination in results_dict.keys():
                OD[o][d] = results_dict[destination]
            # if we can't find a destination node in the results dict, that node cannot be reached from our origin - return the fail value
            else:
                OD[o][d] = fail_value

    # we return our completed, disrupted OD matrix as a numpy array
    return OD

def SummariseOD(OD, fail_value, demand, baseline, GDP_per_capita):
    # this is, by far and away, the most complicated function in this script, and possibly the most complicated I have ever written.
    # this function generates all of the summary statistics we need for the network, based on a passed-in OD matrix.
    # OD - a numpy array of the travel times between all origins and destinations
    # fail_value - the value returned if a trip cannot be completed
    # baseline - the OD matrix when no roads are disrupted. We compare this object to the current OD matrix, which reflects the shortest paths post network disruption
    # demand - this is a numpy matrix of the same dimensions as the OD matrix, which aims to describe the demand between for trips between any node pairs
    # GDP_per_capita - float, GDP per capita

    # set up an empty answers list for population later
    ans = []
    ## calculate trips destroyed by masking the OD matrix for any trips taking longer than the fail value, minus 1.
    # this guarantees that all trips with a time value greater than or equal to the failvalue are chosen.
    # this implies that the SHORTEST path between selected node pairs MUST have used at least one disabled road
    masked_OD = np.ma.masked_greater(OD, value = (fail_value - 1))

    # we do the same thing to the baseline matrix (some trips might not be completable even under perfect circumstances - e.g. island to mainland trips)
    masked_baseline = np.ma.masked_greater(baseline, value = (fail_value - 1))

    # we mask the OD matrix by the baseline mask - to take out the trips that can't be completed even under perfect circumstances
    adj_time = np.ma.masked_array(OD,masked_baseline.mask)

    # we do the same thing to the demand matrix as well for consistency's sake
    masked_demand = np.ma.masked_array(demand, masked_baseline.mask)

    # we calculate the total number of trips undertaken by summing the masked demand matrix
    total_trips = masked_demand.sum()

    # if we use masked OD mask on the demand array, masked for 'bad' trips in perfect circumstances,
    # we can find the trips which only now, post network disruption, can no longer be completed.
    isolated_trips = np.ma.masked_array(masked_demand,~masked_OD.mask)

    # we sum this and slap it into a variable like so
    isolated_trips_sum = isolated_trips.sum()

    # every thing else is a trip that may or may not have been disrupted (but can still be completed).
    # Note we use the INVERSE of the masked_OD mask above, and the straight one here (remember it was masked on the (fail_value - 1))
    # here we are generating the DEMAND matrix first, and then summing it to get total number of 'valid' but potentially disrupted trips.
    potentially_disrupted_trips = np.ma.masked_array(masked_demand,masked_OD.mask)
    potentially_disrupted_trips_sum = potentially_disrupted_trips.sum()

    # Try to calculate the percentage of trips that are isolated. If you can't, then it is by definition 0.
    try:
        pct_isolated = (isolated_trips_sum / total_trips)
    except:
        pct_isolated  = 0

    # we get the corrolary of the array above - the ORIGINAL time of the trips POTENTIALLY disrupted in this scenario.
    potentially_disrupted_trips_original_time = np.ma.masked_array(masked_baseline, masked_OD.mask)

    # we generate a matrix of the actual disruption per the below, called delta-time.
    # This is the time delta incurred by network disruption in this scenario
    delta_time_OD = (masked_OD - potentially_disrupted_trips_original_time)

    # similarly, we can get an average time disruption by summing these matrices and dividing by the original time
    # remember here the 'potentially_disrupted_trips' is the demand matrix.
    average_time_disruption = (delta_time_OD * potentially_disrupted_trips).sum() / potentially_disrupted_trips.sum()

    # we generate a 'fractionalized OD matrix' - by dividing each new time by the original time.
    # Trips now expressed in ratio format (1 = no disruption, >1 = disruption)
    frac_OD = masked_OD / potentially_disrupted_trips_original_time

    # this function identifies the fraction of trips disrupted above a certain time threshold x
    def PctDisrupt(x, frac_OD, demand):
        masked_frac_OD = np.ma.masked_inside(frac_OD, 1, (1+x))
        m_demand = np.ma.masked_array(demand, masked_frac_OD.mask)
        return ((m_demand.sum()) / (demand.sum()))

    # here we run this above disruption function for a value of 30% (x = 0.3)
    pct_thirty_plus = PctDisrupt(0.3, frac_OD, potentially_disrupted_trips)


    def surplus_loss(e, C2, C1, D1):
        # this function works out, at the aggregate economy level, the approximate surplus lost from the network disruption.
        # e - the elasticity of demand w.r.t. trip cost
        # C2 - cost point post disruption
        # C1 - cost point pre disruption
        # D1 - original demand
        # using these four points, we can calculate D2 - the new demand for the
        # trip post disruption - and therefore the size of the surplus loss triangle

        # work out Y intercept of the demand line using known coordinate pair and elasticity
        Y_intercept_max_cost = C1 - (e * D1)

        # work out the Y-axis position - are we at the intercept, or C2?
        C2 = np.minimum(C2, Y_intercept_max_cost)

        # change in trip cost due to disruption is new cost minus original cost
        delta_cost = C2 - C1

        # change in demand is change in cost divided by elasticity
        delta_demand = (delta_cost / e)

        # new demand will be the change in demand plus the original demand (delta_demand often negative)
        D2 = (D1 + delta_demand)

        # use this coordinate to identify surplus loss.
        surplus_loss_ans = ((delta_cost * D2) + ((delta_cost * -delta_demand) / 2))

        # sum up the total surplus loss - this answer an area
        total_surp_loss = surplus_loss_ans.sum()

        # original area of 'surplus' generate by making the trip
        triangle = (D1 * (Y_intercept_max_cost - C1) ) / 2

        # work out the pct. reduction in surplus by dividing one area by the other.
        total_pct_surplus_loss = total_surp_loss / triangle.sum()

        # return the 'absolute' surplus lost, and the pct of original surplus lost.
        return total_surp_loss, total_pct_surplus_loss

    # Here, we convert the OD matrix of adjusted times into a $ figure by using gdp per capita, and assuming time = money
    adj_cost = (adj_time * GDP_per_capita) / (365 * 8 * 3600)

    # we do the same for baseline trip times as well for a point of reference.
    baseline_cost = (masked_baseline * GDP_per_capita) / (365 * 8 * 3600)

    # now, we apply the above function for a low and high bound of the elasticity of demand w.r.t. cost (see J. Rozenberg, LAC DRM for details)
    total_surp_loss_e1, total_pct_surplus_loss_e1 = surplus_loss(-0.15, adj_cost, baseline_cost, masked_demand)
    total_surp_loss_e2, total_pct_surplus_loss_e2 = surplus_loss(-0.36, adj_cost, baseline_cost, masked_demand)

    # end of summarization function - return all scenario summary variables of interest
    return pct_isolated, average_time_disruption, pct_thirty_plus, total_surp_loss_e1, total_pct_surplus_loss_e1, total_surp_loss_e2, total_pct_surplus_loss_e2

def Main(passed_dict):
    ### runs all sub functions - generates altered network, simulates journeys,
    ### and then summarises results. Returns a results dictionary for each scenario run
    # this function is the main workhorse function, performed for every scenario. It is run in parallel by multiple threads in the country_level_criticality function.

    # note that all of the variables required by main are actually bundled into a dictionary called 'passed_dict'.
    # we do this because traditional pythonic multiprocessing does not like the function being parallelized having more than one argument.
    # to get around this, we bundle as many arguments as we like into a dict, pass that as the single argument instead, and then unpack it the other end!

    # unpack our base graph Object
    G = passed_dict['graph']

    # unpack i, the iteration or scenario number
    i = passed_dict['iterno']

    # unpack the fraction of edges to be destroyed
    edge_frac = passed_dict['edge_frac']

    # this is a simple mechanism for keeping track of how many scenarios have been processed by this thread so far.
    if i % 100 == 0:
        print('iteration number {} for edge frac {}'.format(i, edge_frac))

    # break out the list of origin nodes (unchanged from scenario to scenario)
    origins = passed_dict['origins']

    # break out the list of destination nodes (unchanged from scenario to scenario)
    destinations = passed_dict['destinations']

    # break out the fail value
    fail_value = passed_dict['fail_value']

    # break out the demand matrix (unchanged in each scenario - baseline demand)
    demand = passed_dict['demand']

    # unpack the OD travel time matrix, without any adjustments for closed / unpassable roads
    baseline = passed_dict['baseline']

    # unpack the random seed
    seed = passed_dict['seed']

    # unpack the GDP per capita for the country under analysis. Used only for the surplus calculation in SummariseOD
    GDP_per_capita = passed_dict['GDP_per_capita']

    # Perform the random edge destruction on the passed in graph, get back the adjusted network, G_adj
    G_adj, destroy_list = RandomDestruction(G, edge_frac, fail_value, seed)

    # calculate the OD matrix for this scenario, with the adjusted graph, using our origins, destinations, and fail value
    OD = Calculate_OD(G_adj, origins, destinations, fail_value)

    # apply the SummariseOD function on this adjusted OD matrix, with auxillary vars passed in
    pct_isolated, average_time_disruption, pct_thirty_plus, total_surp_loss_e1, total_pct_surplus_loss_e1, total_surp_loss_e2, total_pct_surplus_loss_e2 = SummariseOD(OD, fail_value, demand, baseline, GDP_per_capita)

    # return the results for this scenario expressed as a dictionary, with a few extra
    # bits of info (e.g. i, the scenario number, and the list of edges destroyed)
    return {'iteration':i,
            'nwk_pct_destroyed':edge_frac,
            'pct_journeys_isolated':pct_isolated,
            'average_time_disrupt':average_time_disruption,
            'disrupted_30_pct_plus': pct_thirty_plus,
            'total_surp_loss_e1':total_surp_loss_e1,
            'total_pct_surplus_loss_e1':total_pct_surplus_loss_e1,
            'total_surp_loss_e2':total_surp_loss_e2,
            'total_pct_surplus_loss_e2':total_pct_surplus_loss_e2,
            'edges_destroyed': destroy_list
           }

# this is the main handler function. It is run once per country, and calls the Main function for each scenario.
def country_level_criticality(cou,code,continent,test):

    # print country code to console to ensure we are working on the correct / expected country
    print('current country: %s' % code)

    # here we extract the folder path to the input data required for this analysis from the config.json file
    data_path = load_config()['paths']['data']

    # same as above but the folder path for all of the road networks (which, for C.Fox, were stored separately)
    road_path = os.path.join(load_config()['paths']['road_networks'], code)

    # write path for outputs
    wpath = os.path.join(data_path,'country_networks', code)

    # if a folder at the write path doesn't exist yet, make one.
    # after we've checked that, if a folder called analysis doesn't exist inside the write path folder path, make one
    if not os.path.exists(os.path.join(wpath,'analysis')):
        if not os.path.exists(wpath):
            os.mkdir(wpath)
        os.mkdir(os.path.join(wpath,'analysis'))

    # ### Section 1: Pre Graph Theory:

    # this section could be EASILY modified to use a pickle, and should be.

    # grab the nodes from the folder path
    nodes = pd.read_csv(os.path.join(road_path,'output','{}_processed_nodes.csv'.format(code)))

    # build a geometry field in the GeoDataFrame from the x and y coordinates for each node
    nodes['geometry'] = nodes.apply(lambda z: Point(z.x, z.y), axis=1)

    # make a geodataframe from the nodes
    nodes = gpd.GeoDataFrame(nodes)

    # load edges from road_path
    edges = pd.read_csv(os.path.join(road_path,'output','{}_processed_edges.csv'.format(code)))

    # ensure we have a field for each edge called 'highway'
    edges.rename({'infra_type':'highway'},axis=1,inplace=True)

    # make the geometry field equal to the loaded WKT shapely geometry currently stored as string in the 'geometry' field
    s = edges['geometry']
    l = s.apply(shapely.wkt.loads)
    edges.geometry=l

    # with these bit assembled, make the edges GeoDataFrame
    edges = gpd.GeoDataFrame(edges)

    # send the edges to a shapefile for later analysis / checking
    edges.to_file(os.path.join(wpath,'analysis','{}_simplev2_edges.shp'.format(code)))

    # ### Section 2: Import Shapefile as Graph Object, prepare for analysis

    # the shape of the graph needs to be fairly specific for this analysis. Several edge properties need to be properly in place:
    # - there needs to be an 'edge_id' fields that is an integer, counting upwards
    # - each edge needs a 'time' field, measured in seconds
    # - for every edge linking node u to v, there needs to be an edge linking v to u.
    #   Networkx's default import of a shapefile imports a unidirectional DiGraph - so be sure to make sure you can get from v to u as well as u to v
    # - the default import labels the nodes with their geometry position. create a new integer mapping for the node IDs. This prevents painful object type confusion later in the analysis

    # Here we read in the edges shapefile as a networkx object (WARNING - in future, read a gpickle - much cleaner)
    graph = nx.read_shp(os.path.join(wpath,'analysis','{}_simplev2_edges.shp'.format(code)))

    # here we label up the nodes
    z = 0
    mapping = {}
    for u, data in graph.nodes(data = True):
        data['x'] = u[0]
        data['y'] = u[1]
        mapping[u] = z
        z += 1

    # Relabel nodes with integer value
    graph = nx.relabel_nodes(graph, mapping, copy=True)

    # Convert to MultiDigraph - add reflected edges
    new_edge_bucket = []

    # set up the integer edge ID counter field.
    edge_id = 0
    for u,v, data in graph.edges(data = True):
        data['edge_id'] = edge_id
        # here, I am specifying the inverse of the edge currently under consideration
        # and adding that with the same 'id' property as the current edge. So, the reflected edge will have the same ID. important!
        new_edge = (v,u, data)

        # add the new edge to the new edge bucket
        new_edge_bucket.append(new_edge)
        edge_id += 1

    # add everything in the new edge bucket to the graph
    graph.add_edges_from(new_edge_bucket)

    # set up an edge GDF for future reference later
    graph_df = gn.edge_gdf_from_graph(graph)

    # Add time attribute to edges. All speeds in kmph
    speed_d = {
    'motorway':80,
    'motorway_link': 65,
    'trunk': 60,
    'trunk_link':50,
    'primary': 50, # kmph
    'primary_link':40,
    'secondary': 40, # kmph
    'secondary_link':30,
    'tertiary':30,
    'tertiary_link': 20,
    'unclassified':20,
    'residential': 20,  # kmph
    }

    # use a GOSTNets function to add on the time factor. 'highway' is the default edge property used to assign speeds and times
    graph_time = gn.convert_network_to_time(graph, distance_tag = 'length', graph_type = 'drive', speed_dict = speed_d, factor = 1000)

    # Section 3: Generate O-D points

    # Import the boundaries and population files

    # bring in the boundaries file - global - to be subsetted later
    bounds = os.path.join(data_path,'input_data','global_regions_v2.shp')

    #set global file paths for worldpop. continent variable decides which is used
    if (continent == 'Central-America') or (continent == 'South-America') :
        wp = os.path.join(data_path,'Worldpop','LAC_PPP_2015_adj_v2.tif')
    elif (continent == 'Africa'):
        wp = os.path.join(data_path,'Worldpop','AFR_PPP_2015_adj_v2.tif')
    elif (continent == 'Asia'):
        wp = os.path.join(data_path,'Worldpop','Asia_PPP_2015_adj_v2.tif')
    elif (continent == 'Europe'):
        wp = os.path.join(data_path,'Worldpop','EUROPOP_WGS84.tif')

    # due to a few islands not included in the global worldpop data, we need to
    # extract their data individually
    islands = ['FJI','KIR','MHL','FSM','PLW','PNG','WSM','SLB','TON','VUT']

    # is we are looking at an insland, use the unzip_worldpop function to grab the WP info
    if code in islands:
        temp_path = os.path.join(data_path,'Worldpop','temp_{}'.format(code))

        # make a temporary folder
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

        # unzip the world pop into the temp folder
        unzip_worldpop(code,data_path,temp_path)
        wp = os.path.join(temp_path,'popmap15adj.tif')
        if code == 'KIR':
            wp = os.path.join(temp_path,'popmap15adj_lzw.tif')
        elif code == 'TON':
            wp = os.path.join(temp_path,'popmap15.tif')
        wp_out = os.path.join(temp_path,'popmap15adj_wgs84.tif')
        os.system('gdalwarp -t_srs EPSG:4326 -tr 0.018 0.018 -overwrite '+wp+' '+wp_out)
        wp = wp_out

    # special instructions for Mexico and Ukraine - specific to World Pop file formatting
    if code == 'MEX':
        wp = os.path.join(data_path,'Worldpop','MEX_ppp_v2c_2015_UNadj.tif')
    elif code == 'UKR':
        wp = os.path.join(data_path,'Worldpop','UKR_ppp_v2b_2015_UNadj.tif')

    # read the global boundary file, use the passed in ISO code to pick out the correct shape.
    # note that this is at the GADM 2 level - i.e. regions / districts. Multiple polygons per ISO code.
    b = gpd.read_file(bounds)
    b=b[b['GID_0']==code]

    # set ID and area values from the .loc'd geodataframe
    b['ID'] = b.index
    b['area'] = b.area

    # send a copy of this to file for later reference.
    b.to_file(os.path.join(wpath,'analysis','{}.shp'.format(code)))

    country_bounds = b

    # here, we use the zonal stats function against the world pop path, wp, to summarize
    # the worldpop in each admin region inside the country of interest
    b['sum'] = zonal_stats(country_bounds, wp, stats="sum")

    # we sum the sums to get the total population sum of the country across all polygons
    b['sum'] = b['sum'].apply(lambda x: x['sum'])

    # change the geometry to the centroid of each admin region
    b['geometry'] = b.centroid

    # TP is the total number of points we want for both origins and destinations
    TP = 100

    # the analysis is set to choose 50% of the admin 2 regions based on population (top50),
    # and 50% of the points based on the area of the admin regions (excl. the Top 50 most populous regions)
    pop_share  = 0.5
    size_share = 0.5

    # here, we effect the above logic against c, a copy of b, the admin regions GeoDataFrame
    c = b.copy()
    c = c.sort_values(by = 'sum', ascending = False)
    d = c.copy()

    # take total points X pop_share of the dataframe using a row-slice
    c = c[:int(TP*pop_share)]

    # define d as the 'other' admin regions not in the above frames
    d = d[len(c):]

    # sort these by area
    d = d.sort_values(by = 'area', ascending = False)

    # select total points X area_share of these for origins too
    d = d[:int(TP*size_share)]

    # make a new frame, a copy of the original
    e = b.copy()

    # slice it for ONLY those points in either c or d.
    e = e.loc[(b.ID.isin(d.ID)) | (b.ID.isin(c.ID))]

    # match these points in e onto the graph, returning the nearest node to each point.
    e = gn.pandana_snap(graph_time, e)

    # here we copy and subset the smaller dataset to just the nearest node and the sum of the population in that polygon
    f = e.copy()
    f = f[['sum','NN']]

    # here i do some shenanigans to get the unique origins (in terms of nodes, NOT centroids) and
    # associated snapped population into list objects. Don't look to closely, it works
    f['dummy'] = 1
    f = f.set_index(['NN','dummy'])
    f = f.sum(level = 'NN')
    unique_origins = list(f.index)
    pop = list(f['sum'])

    # Here, we build the demand matrix to accompany the baseline OD matrix.
    # We normalize the number of trips to 'maxtrips'
    maxtrips = 100

    # this parameter allows you to tweak the aggression of the distance decay behaviour on the demand between any two node pairs
    dist_decay = 1

    # here we define our fail value, which is returned if a journey cannot be completed between A and B
    # (implies more than one subgraph or road disrupted)
    fail_value = 999999999999

    # we set up a blank numpy array for the demand, same length as the origins and destinations
    demand = np.zeros((len(unique_origins), len(unique_origins)))

    # we calculate the baseline OD matrix, which has in it the travel time between node pairs assuming no disruption
    shortest_time = Calculate_OD(graph_time, unique_origins, unique_origins, fail_value)

    # Calculate demand between each origin and destination
    for o in range(0, len(unique_origins)):
        for d in range(0, len(unique_origins)):
            if o == d:
                # do not insert demand down the spine - no trips where origin = destination
                demand[o][d] = 0
            else:
                # normalize the current travel time versus the largest travel time between nodes in the matrix
                normalized_dist = shortest_time[o][d] / shortest_time.max()

                # here, demand is a function of the product of the population of the origin and
                #  the destination - but reduced by the distance between them. 'Gravity demand'
                demand[o][d] = ((pop[o] * pop[d]) * np.exp(-1 * dist_decay * normalized_dist))

    # we noralize the matrix to the number of maxtrips
    demand = ((demand / demand.max()) * maxtrips)

    # we round up - to ensure each journey is made at least once
    demand = np.ceil(demand).astype(int)


    # CONTROLS
    # here, we define the number of random scenarios tested for each level of disruption set in edge_fracs
    iterations = 500

    # This is the % of the network that will be destroyed in a given scenario
    edge_fracs = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]

    # set origins and destinations equal to the unique listed of snapped origins and destinations
    origins = unique_origins
    destinations = unique_origins

    # open a blank results list
    results = []
    baseline = shortest_time

    # identify trips that already can't be completed, absent network disruptions
    masked_baseline = np.ma.masked_greater(baseline, value = (fail_value - 1))

    # Here is set up the number of threads, which is equal to the CPU count.
    # You can limit this to a lower (but not higher number) than the CPU count
    threads = multiprocessing.cpu_count()

    # Read in the GDP dataframe
    GDP_df = pd.read_csv(os.path.join(data_path, 'input_data','gdp.csv'))

    # Pick out the GDP per capita for the current country
    GDP_per_capita = GDP_df.latest_available.loc[GDP_df['Country Code'] == code].iloc[0]
    print('GDP per capita:', GDP_per_capita)

    # set up some timing
    start = time.time()
    print("Commencing Main with {} thread(s) for {}. Start time: {}" .format(threads, cou, time.ctime()))

    # in this block we manufacture the dictionary of variables that is passed to each thread when processing a scenario.
    # we load all these many dictionaries into a massive list, which we then distribute amongst the avaialble threads
    d = []

    # we generate a different random scenario for each scenario by iterating up
    # through the random seed for every dictionary of parameters that we generate
    seed = 1
    # for each disruption fraction...
    for edge_frac in edge_fracs:
        # generate 'iterations' iterations...
        for i in range(0,iterations):
            # build the dict of necessary variables..
              l = {'graph':graph_time,
                  'iterno':i,
                  'edge_frac':edge_frac,
                  'origins':origins,
                  'destinations':destinations,
                  'fail_value':fail_value,
                  'demand':demand,
                  'baseline':baseline,
                  'seed':seed,
                  'GDP_per_capita':GDP_per_capita}
              # append this to the monster list
              d.append(l)
              seed+=1

    # set up a multiprocessing Pool
    with Pool(threads) as pool:

        # create an object called results which is the aggregate of all of the
        # threads and their resoluts when they execute the function 'Main'
        results = pool.map(Main,d,chunksize=1)

    # timing
    elapsed = (time.time() - start)
    print('finished {} after {} seconds, {} scenarios, or c. {} seconds per scenario'.format(cou,elapsed, iterations*len(edge_fracs), (elapsed / (iterations*len(edge_fracs)))))

    # make a dataframe of the results
    l = pd.DataFrame(results)
    l = l[['nwk_pct_destroyed',
    'iteration',
    'pct_journeys_isolated',
    'average_time_disrupt',
    'disrupted_30_pct_plus',
    'total_surp_loss_e1',
    'total_pct_surplus_loss_e1',
    'total_surp_loss_e2',
    'total_pct_surplus_loss_e2',
    'edges_destroyed']]

    # create the outpath, write file to disk. Done!
    outpath = os.path.join(data_path,'country_networks',code,'analysis')
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    l.to_csv(os.path.join(os.path.join(outpath,'{}_criticality_res.csv'.format(cou))))

# Pythonically, this is where the script starts when called from the command line.
# the primary job of main is to launch the function country_level_criticality
if __name__ == '__main__':

    # load the datapath...again...
    data_path = load_config()['paths']['data']

    # load the global shapefile
    global_countries = gpd.read_file(os.path.join(data_path,'input_data','global_countries.shp'))

    # make a list of the continents, iso codes and countries
    iso3 = list(global_countries.ISO_3digit)
    countries = list(global_countries.NAME_0)
    continent = list(global_countries.continent)

    # set up testing - if testing, turn test to 1. It will just do the
    # first country by subsetting the list inputs appropriately
    test = 0
    if test == 1:
        iso3 = iso3[:1]
        countries = countries[:1]
        continent = continent[:1]

    # this executes the country_level_criticality function for every ISO code:
    for code,cou,continent in zip(iso3,countries,continent):
        country_level_criticality(cou,code,continent,test)
