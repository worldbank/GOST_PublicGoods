
# coding: utf-8
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
sys.path.append(os.getcwd())
from utils import load_config
import multiprocessing

sys.path.append(r'/home/wb493355/data/criticality II/scripts/')
import GOSTnet as gn

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
    ### function randomly sets the traverse time on chose edges to the fail_value, simulating road unavailability
    ### proportion of total edges disrupted is the edge_frac

    edgeid = []
    for u,v, data in G.edges(data = True):
        edgeid.append(data['edge_id'])
    edgeid = list(set(edgeid))

    num_to_destroy = math.floor(len(edgeid) * (edge_frac / 100))
    np.random.seed(seed = seed)
    np.random.shuffle(edgeid)
    destroy_list = edgeid[:num_to_destroy]

    G_adj = G.copy()
    destroy_count = 0
    for u, v, data in G_adj.edges(data = True):
        if data['edge_id'] in destroy_list:
            data['time'] = fail_value
            destroy_count +=1
    return G_adj, destroy_list

def Calculate_OD(G, origins, destinations, fail_value, weight = 'time'):
    ### returns numpy O-D array for shortest travel time between each origins and destination.
    ### travel time between o and d is 'answers[o][d]'
    ### incompletable trips return the fail_value

    OD = np.zeros((len(origins), len(destinations)))

    for o in range(0, len(origins)):
        origin = origins[o]
        results_dict = nx.single_source_dijkstra_path_length(G, origin, cutoff = None, weight = weight)

        for d in range(0, len(destinations)):
            destination = destinations[d]
            if destination in results_dict.keys():
                OD[o][d] = results_dict[destination]
            else:
                OD[o][d] = fail_value

    return OD

def SummariseOD(OD, fail_value, demand, baseline, GDP_per_capita):

    ans = []

    ### Function returns the % of total trips between origins and destinations that exceed fail value

    ## calculate trips destroyed
    masked_OD = np.ma.masked_greater(OD, value = (fail_value - 1))
    masked_baseline = np.ma.masked_greater(baseline, value = (fail_value - 1))
    adj_time = np.ma.masked_array(OD,masked_baseline.mask)
    masked_demand = np.ma.masked_array(demand, masked_baseline.mask)

    total_trips = masked_demand.sum()

    isolated_trips = np.ma.masked_array(masked_demand,~masked_OD.mask)
    isolated_trips_sum = isolated_trips.sum()

    potentially_disrupted_trips = np.ma.masked_array(masked_demand,masked_OD.mask)
    potentially_disrupted_trips_sum = potentially_disrupted_trips.sum()

    try:
        pct_isolated = (isolated_trips_sum / total_trips)
    except:
        pct_isolated  = 0

    ## Set up for ratio-based calculations

    fail_ratio = 50.0 #((fail_value-1) / baseline.max()) # headlimit above which trip destoryed

    potentially_disrupted_trips_original_time = np.ma.masked_array(masked_baseline, masked_OD.mask)
    delta_time_OD = (masked_OD - potentially_disrupted_trips_original_time)
    average_time_disruption = (delta_time_OD * potentially_disrupted_trips).sum() / potentially_disrupted_trips.sum()

    frac_OD = masked_OD / potentially_disrupted_trips_original_time

    def PctDisrupt(x, frac_OD, demand):
        masked_frac_OD = np.ma.masked_inside(frac_OD, 1, (1+x))
        m_demand = np.ma.masked_array(demand, masked_frac_OD.mask)
        return ((m_demand.sum()) / (demand.sum()))

    pct_thirty_plus = PctDisrupt(0.3, frac_OD, potentially_disrupted_trips)

    # Flexing demand with trip cost
    def surplus_loss(e, C2, C1, D1):

        Y_intercept_max_cost = C1 - (e * D1)

        C2 = np.minimum(C2, Y_intercept_max_cost)

        delta_cost = C2 - C1

        delta_demand = (delta_cost / e)

        D2 = (D1 + delta_demand)

        surplus_loss_ans = ((delta_cost * D2) + ((delta_cost * -delta_demand) / 2))

        triangle = (D1 * (Y_intercept_max_cost - C1) ) / 2

        total_surp_loss = surplus_loss_ans.sum()

        total_pct_surplus_loss = total_surp_loss / triangle.sum()

        return total_surp_loss, total_pct_surplus_loss

    adj_cost = (adj_time * GDP_per_capita) / (365 * 8 * 3600)
    baseline_cost = (masked_baseline * GDP_per_capita) / (365 * 8 * 3600)

    total_surp_loss_e1, total_pct_surplus_loss_e1 = surplus_loss(-0.15, adj_cost, baseline_cost, masked_demand)
    total_surp_loss_e2, total_pct_surplus_loss_e2 = surplus_loss(-0.36, adj_cost, baseline_cost, masked_demand)

    return pct_isolated, average_time_disruption, pct_thirty_plus, total_surp_loss_e1, total_pct_surplus_loss_e1, total_surp_loss_e2, total_pct_surplus_loss_e2

    """
    def pct_calc(frac_OD, demand, pct, fail_ratio):
        # n.b. - only actually disrupted trips. trips with a ratio of 1 remain the same
        masked_OD = np.ma.masked_outside(frac_OD, 1.001, float((1 + pct)))
        disrupted = (np.ma.masked_array(masked_demand, masked_OD.mask))

        try:
            pct_disrupt = (disrupted.sum() / potentially_disrupted_trips_sum) * 100
        except:
            pct_disrupt  = 0

        return pct_disrupt
    """

def Main(passed_dict):
    ### runs all sub functions - generates altered network, simulates journeys,
    ### and then summarises results. Returns a results dictionary for each scenario run

    G = passed_dict['graph']
    i = passed_dict['iterno']
    edge_frac = passed_dict['edge_frac']
    if i % 100 == 0:
        print('iteration number {} for edge frac {}'.format(i, edge_frac))
    origins = passed_dict['origins']
    destinations = passed_dict['destinations']
    fail_value = passed_dict['fail_value']
    demand = passed_dict['demand']
    baseline = passed_dict['baseline']
    seed = passed_dict['seed']
    GDP_per_capita = passed_dict['GDP_per_capita']

    G_adj, destroy_list = RandomDestruction(G, edge_frac, fail_value, seed)

    OD = Calculate_OD(G_adj, origins, destinations, fail_value)

    pct_isolated, average_time_disruption, pct_thirty_plus, total_surp_loss_e1, total_pct_surplus_loss_e1, total_surp_loss_e2, total_pct_surplus_loss_e2 = SummariseOD(OD, fail_value, demand, baseline, GDP_per_capita)

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

# ### Section 6: Run Main Function

def country_level_criticality(cou,code,continent,test):
#    try:
    print('current country: %s' % code)
    data_path = load_config()['paths']['data']
    road_path = os.path.join(load_config()['paths']['road_networks'], code)
    wpath = os.path.join(data_path,'country_networks', code)
    if not os.path.exists(os.path.join(wpath,'analysis')):
        if not os.path.exists(wpath):
            os.mkdir(wpath)
        os.mkdir(os.path.join(wpath,'analysis'))

    # ### Section 1: Pre Graph Theory:
    nodes=pd.read_csv(os.path.join(road_path,'output','{}_processed_nodes.csv'.format(code)))
    nodes['geometry'] = nodes.apply(lambda z: Point(z.x, z.y), axis=1)
    nodes = gpd.GeoDataFrame(nodes)
    edges=pd.read_csv(os.path.join(road_path,'output','{}_processed_edges.csv'.format(code)))
    edges.rename({'infra_type':'highway'},axis=1,inplace=True)
    s=edges['geometry']
    l=s.apply(shapely.wkt.loads)
    edges.geometry=l
    edges = gpd.GeoDataFrame(edges)
    #nodes.to_file(os.path.join(wpath,'analysis','{}_simplev2_nodes.shp'.format(code)))
    edges.to_file(os.path.join(wpath,'analysis','{}_simplev2_edges.shp'.format(code)))

    # ### Section 2: Import Shapefile as Graph Object, prepare for analysis

    # the shape of the graph needs to be fairly specific for this analysis. Several edge properties need to be properly in place:
    # - there needs to be an 'edge_id' fields that is an integer, counting upwards
    # - each edge needs a 'time' field, measured in seconds
    # - for every edge linking node u to v, there needs to be an edge linking v to u. Networkx's default import of a shapefile imports a unidirectional DiGraph - so be sure to make sure you can get from v to u as well as u to v
    # - the default import labels the nodes with their geometry position. create a new integer mapping for the node IDs. This prevents painful object type confusion later in the analysis

    graph = nx.read_shp(os.path.join(wpath,'analysis','{}_simplev2_edges.shp'.format(code)))
    z = 0
    mapping = {}
    for u, data in graph.nodes(data = True):
        data['x'] = u[0]
        data['y'] = u[1]
        mapping[u] = z
        z += 1

    # #### Relabel nodes with integer value

    graph = nx.relabel_nodes(graph, mapping, copy=True)

    # #### Convert to MultiDigraph - add reflected edges

    new_edge_bucket = []
    edge_id = 0
    for u,v, data in graph.edges(data = True):
        data['edge_id'] = edge_id
        new_edge = (v,u, data)
        new_edge_bucket.append(new_edge)
        edge_id += 1
    graph.add_edges_from(new_edge_bucket)

    graph_df = gn.edge_gdf_from_graph(graph)

    # #### Add time attribute to edges

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

    graph_time = gn.convert_network_to_time(graph, distance_tag = 'length', graph_type = 'drive', speed_dict = speed_d, factor = 1000)

    ## TO DO - add in Julie Funcs here ##

# ### Section 3: Generate O-D points

# #### Import the boundaries and population files

    # ### Section 3: Generate O-D points

    # #### Import the boundaries and population files

    bounds = os.path.join(data_path,'input_data','global_regions_v2.shp')

        #set global file paths for worldpop
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
    if code in islands:
        temp_path = os.path.join(data_path,'Worldpop','temp_{}'.format(code))

        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

        unzip_worldpop(code,data_path,temp_path)
        wp = os.path.join(temp_path,'popmap15adj.tif')
        if code == 'KIR':
            wp = os.path.join(temp_path,'popmap15adj_lzw.tif')
        elif code == 'TON':
            wp = os.path.join(temp_path,'popmap15.tif')
        wp_out = os.path.join(temp_path,'popmap15adj_wgs84.tif')
        os.system('gdalwarp -t_srs EPSG:4326 -tr 0.018 0.018 -overwrite '+wp+' '+wp_out)
        wp = wp_out

    if code == 'MEX':
        wp = os.path.join(data_path,'Worldpop','MEX_ppp_v2c_2015_UNadj.tif')
    elif code == 'UKR':
        wp = os.path.join(data_path,'Worldpop','UKR_ppp_v2b_2015_UNadj.tif')

    b = gpd.read_file(bounds)
    b=b[b['GID_0']==code]
    b['ID'] = b.index
    b['area'] = b.area
    b.to_file(os.path.join(wpath,'analysis','{}.shp'.format(code)))
    country_bounds = os.path.join(wpath,'analysis','{}.shp'.format(code))

    b['sum'] = zonal_stats(country_bounds, wp, stats="sum")
    b['sum'] = b['sum'].apply(lambda x: x['sum'])
    #b.to_file(os.path.join(os.path.join(wpath,'analysis','{}_poly_pop.shp'.format(cou))))
    b['geometry'] = b.centroid

    TP = 100
    pop_share  = 0.5
    size_share = 0.5

    c = b.copy()
    c = c.sort_values(by = 'sum', ascending = False)
    d = c.copy()
    c = c[:int(TP*pop_share)]
    d = d[len(c):]
    d = d.sort_values(by = 'area', ascending = False)
    d = d[:int(TP*size_share)]
    e = b.copy()
    e = e.loc[(b.ID.isin(d.ID)) | (b.ID.isin(c.ID))]
    #e.to_file(os.path.join(os.path.join(wpath,'analysis','{}_top.shp'.format(cou))))

    e = gn.pandana_snap(graph_time, e)

    f = e.copy()
    f = f[['sum','NN']]
    f['dummy'] = 1
    f = f.set_index(['NN','dummy'])
    f = f.sum(level = 'NN')
    unique_origins = list(f.index)
    pop = list(f['sum'])
    maxtrips = 100
    dist_decay = 1

    fail_value = 999999999999

    demand = np.zeros((len(unique_origins), len(unique_origins)))

    shortest_time = Calculate_OD(graph_time, unique_origins, unique_origins, fail_value)

    for o in range(0, len(unique_origins)):
        for d in range(0, len(unique_origins)):
            if o == d:
                demand[o][d] = 0
            else:
                normalized_dist = shortest_time[o][d] / shortest_time.max()
                demand[o][d] = ((pop[o] * pop[d]) * np.exp(-1 * dist_decay * normalized_dist))

    demand = ((demand / demand.max()) * maxtrips)
    demand = np.ceil(demand).astype(int)

    iterations = 500
    edge_fracs = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
    origins = unique_origins
    destinations = unique_origins

    results = []
    baseline = shortest_time
    masked_baseline = np.ma.masked_greater(baseline, value = (fail_value - 1))

    threads = multiprocessing.cpu_count()

    GDP_df = pd.read_csv(os.path.join(data_path, 'input_data','gdp.csv'))
    GDP_per_capita = GDP_df.latest_available.loc[GDP_df['Country Code'] == code].iloc[0]
    print('GDP per capita:', GDP_per_capita)
    start = time.time()
    print("Commencing Main with {} thread(s) for {}. Start time: {}" .format(threads, cou, time.ctime()))

    d = []
    seed = 1
    for edge_frac in edge_fracs:
        for i in range(0,iterations):
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
              d.append(l)
              seed+=1

    with Pool(threads) as pool:
        results = pool.map(Main,d,chunksize=1)
    elapsed = (time.time() - start)

    print('finished {} after {} seconds, {} scenarios, or c. {} seconds per scenario'.format(cou,elapsed, iterations*len(edge_fracs), (elapsed / (iterations*len(edge_fracs)))))

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

    outpath = os.path.join(data_path,'country_networks',code,'analysis')
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    l.to_csv(os.path.join(os.path.join(outpath,'{}_criticality_res.csv'.format(cou))))

if __name__ == '__main__':
    data_path = load_config()['paths']['data']
    global_countries = gpd.read_file(os.path.join(data_path,'input_data','global_countries.shp'))
    #global_countries = global_countries.loc[global_countries.continent == 'Africa'][::-1]

    iso3 = list(global_countries.ISO_3digit)
    countries = list(global_countries.NAME_0)
    continent = list(global_countries.continent)
    test = 0
    if test == 1:
        iso3 = iso3[:1]
        countries = countries[:1]
        continent = continent[:1]
    for code,cou,continent in zip(iso3,countries,continent):
        if code in ['MDG']:
            country_level_criticality(cou,code,continent,test)
