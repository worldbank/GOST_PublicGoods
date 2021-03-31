#!/usr/bin/env python
# coding: utf-8

# # Step 3: Calculate OD

# In[2]:


import os, sys, time, importlib
import osmnx

import geopandas as gpd
import pandas as pd
import networkx as nx
import numpy as np
sys.path.append("../../../GOSTnets")
import GOSTnets as gn

from shapely.geometry import LineString, Point


# In[3]:


# This is a Jupyter Notebook extension which reloads all of the modules whenever you run the code
# This is optional but good if you are modifying and testing source code
#get_ipython().run_line_magic('load_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')


# In[4]:


from GOSTnets.load_traffic2 import *


# In[5]:


# read graph
G = nx.read_gpickle('./sri_lanka_unclean2_w_time_largest_20200616_traffic_mean_speeds.pickle')


# In[6]:


#len(G.edges)


# In[7]:


#gn.example_edge(G, 5)


# ## load origins

# In[8]:


origins = gpd.read_file('./origins_destinations/intl_airport_updated.shp')


# In[9]:


origins


# ## load destinations

# In[10]:


destinations = gpd.read_file('./origins_destinations/cities_top10.shp')


# In[11]:


destinations


# In[12]:


origins_gdf = gn.pandana_snap_c(G, origins, source_crs = 'epsg:4326', target_crs = 'epsg:32644')


# In[13]:


#origins_gdf


# In[14]:


origins_list = list(set(origins_gdf.NN))


# In[15]:


destinations_gdf = gn.pandana_snap_c(G, destinations, source_crs = 'epsg:4326', target_crs = 'epsg:32644')


# In[16]:


#destinations_gdf


# In[17]:


destinations_list = list(set(destinations_gdf.NN))


# In[18]:


#destinations_list


# ## Calculate OD

# In[19]:


# It will use the default weight of 'time'
import time

start = time.time()

OD = gn.calculate_OD(G, origins_list, destinations_list, fail_value = 9999999)

end = time.time()
print(end - start)


# In[20]:


#OD


# In[21]:


OD_df = pd.DataFrame(OD, index = origins_list, columns = destinations_list)


# In[22]:


OD_df


# In[23]:


#OD_df.min(axis=0)


# ### Now we need to find the nearest city for each tourist point, and based on this assign a group of tourist points to each of the 10 cities

# In[24]:


# test
# first column is the destination (tourist point), and the second column is the nearest city
OD_df.idxmin(axis=0)


# In[25]:


#OD_df.idxmin(axis=0).to_csv('./nodes_and_associated_nearest_city.csv')


# In[26]:


#OD_df.idxmin(axis=0).to_frame(0)


# ### takes the min index value of each column, then groups by city (first index (0)) and takes the first entry

# In[27]:


groupby_obj = OD_df.idxmin(axis=0).to_frame(0).groupby(0)[0]


# In[28]:


groupby_obj


# In[29]:


# visualize groupby_obj
#groupby_obj.apply(list)


# In[30]:


#type(groupby_obj.apply(list))


# In[31]:


# a nice way to visualize the groupby_obj
#groupby_obj.describe()


# ### create a dictionary that associates assigned tourist points with each city

# In[32]:


city_tourist_pt_dict = {}
for name, group in groupby_obj:
    #print(group)
    for items in group.iteritems(): 
        #print(items[1])
        if items[1] not in city_tourist_pt_dict:
            city_tourist_pt_dict[items[1]] = [items[0]]
        else:
            #append value to list in dict value
            city_tourist_pt_dict[items[1]].append(items[0])
    #print(type(group))
    #print(group.head(1))

    #print(name)
    #print(city_tourist_pt_dict[group])


# In[33]:


city_tourist_pt_dict


# ## Loop through dictionary in order to do a calculate_OD for each airport's nearest cities

# In[34]:


OD = {}
OD_df_dict = {}
for city,tourist_pts in city_tourist_pt_dict.items():
    OD[city] = gn.calculate_OD(G, [city], tourist_pts, fail_value = 9999999)
    OD_df_dict[city] = pd.DataFrame(OD[city], index = [city], columns = tourist_pts)


# In[35]:


#sample_df = OD_df_dict[4784044133]
#sample_df


# In[36]:


OD_df_dict


# ### time calculateOD time for airport city and its set of city points

# In[37]:


#start = time.time()
#OD[4784044133] = gn.calculate_OD(G, [4784044133], tourist_pts, fail_value = 9999999)
#print(time.time() - start)


# In[38]:


# sum of all shortest routes
#OD[4784044133][0].sum()


# ## Now work on generating routes and visualizing them

# In[39]:


from shapely.ops import linemerge
from itertools import islice


# In[40]:


maxspeed_mean_speeds = {
    'secondary': 50,
    'secondary_link': 45,
    'tertiary': 40,
    'tertiary_link': 40,
    'residential': 25,
    'unclassified': 25,
}

mapbox_mean_speeds = {
    'secondary': 34,
    'secondary_link': 9,
    'tertiary': 25,
    'tertiary_link': 13,
    'residential': 20,
    'unclassified': 20,
}

# In[41]:


def tabulate_edges(route):
    edge_table = []
    route_geometry = LineString()
    for idx in range(0, len(route) - 1):
        #edge_table.append([route[idx], route[idx+1]])
        # look up line
        #print('to node')
        #print(route[idx])
        #print('from node')
        #print(route[idx+1])
        edge_geometry = G.get_edge_data(route[idx],route[idx+1])[0]['geometry']
        # get edge speed
        edge_speed = G.get_edge_data(route[idx],route[idx+1])[0]['speed']
        #print('print edge_speed')
        #print(edge_speed)
        
        # compare edge speed to median speed
        rural_roads_list = ['residential','secondary','secondary_link','tertiary','tertiary_link','unclassified']
        
        edge_infra_type = G.get_edge_data(route[idx],route[idx+1])[0]['infra_type']
        
        edge_length = G.get_edge_data(route[idx],route[idx+1])[0]['length']
        edge_time = G.get_edge_data(route[idx],route[idx+1])[0]['time']
        
        mean_speed = G.get_edge_data(route[idx],route[idx+1])[0]['mean_speed']
        
        try:
            edge_imp_cost = G.get_edge_data(route[idx],route[idx+1])[0]['imp_cost']
        except:
            edge_imp_cost = 0
            pass
        
        if mean_speed > 0:      
            if edge_infra_type in rural_roads_list:
                print('print edge attributes')
                print(G.get_edge_data(route[idx],route[idx+1])[0])
                #print('print improved speed')
                #print(mapbox_traffic_mean_speeds[edge_infra_type])
                #assums that current edge lenght is in km
                new_time_s = (edge_length / maxspeed_mean_speeds[edge_infra_type]) * 3600
                #new_time_s = ((edge_length / 1000) / maxspeed_mean_speeds[edge_infra_type]) * 3600
                #print("print new_time_s")
                #print(new_time_s)
                #print('print old time')
                #print(edge_time)
                edge_savings = edge_time - new_time_s
                # assign savings time
                edge_table.append([route[idx], route[idx+1], edge_savings, edge_imp_cost, edge_length, edge_infra_type, mean_speed, edge_geometry])
                # data['sec_saved'] = data['length'] / new_time
                # data['improvement_cost'] = 174861 * data['length'] / 1000
            else:
                # data['sec_saved'] = 0
                # very high number
                # data['improvement_cost'] = 10000000000
                #print("not a rural road")
                edge_table.append([route[idx], route[idx+1], 0, 0, edge_length, edge_infra_type, mean_speed, edge_geometry])
        else:
            edge_table.append([route[idx], route[idx+1], 0, 0, edge_length, edge_infra_type, mean_speed, edge_geometry])
            
        route_geometry = route_geometry.union(edge_geometry)
        
    #print('print route_geometry')
    #print(route_geometry)
    
    return(edge_table, route_geometry)


# In[42]:


def generate_complete_edges_and_routes(input_df):

    LIMIT = 1000000

    complete_edges = []
    complete_routes = []

    count = 0

    # for origin, row in sample_df.iterrows(): 
    for origin, row in islice(input_df.iterrows(), LIMIT):    
        for destination, value in islice(row.items(), LIMIT):
            try:
                origin = int(origin)
            except:
                pass
            try:
                destination = int(destination)
            except:
                pass

            count = count + 1

            route = nx.dijkstra_path(G, origin, destination, weight = 'time')
            #path_edges = zip(route,route[1:])
            #print('print path_edges')
            #print(list(path_edges))
            # print(route)
            edge_table, route_geometry = tabulate_edges(route)
            #print('print edge_table:')
            #print(edge_table)
            complete_edges = complete_edges + edge_table
            #print('route_time')
            #print(value)
            complete_routes.append([edge_table[0][0], edge_table[-1][1], value, route_geometry])
            #print('edge_table[:-1]')
            #print(edge_table[-1][1])
            
    # convert complete_edges to gdf
    complete = pd.DataFrame(complete_edges, columns = ['o', 'd', 'sec_saved', 'imp_cost', 'length', 'infra_type', 'mean_speed', 'geometry'])
    complete['w'] = 1
    complete_count = complete.groupby(['o','d']).agg(
        {
            'w':"count",
            'sec_saved': 'first',
            'imp_cost': 'first',
            'mean_speed': 'first',
            'length':'first',
            'infra_type':'first',
            'geometry':'first'
        }
    )
    complete_count.reset_index(inplace = True)
    complete_count['o'] = complete_count['o'].astype(str)
    complete_count['d'] = complete_count['d'].astype(str)
    complete_count['weighted_sec_saved'] = complete_count.w * complete_count.sec_saved
    complete_count.sort_values(by=['weighted_sec_saved'], ascending=False)
    complete_count_gdf = gpd.GeoDataFrame(complete_count, crs = 'epsg:4326')
    
    # convert complete_routes to gdf
    complete_routes_df = pd.DataFrame(complete_routes, columns = ['origin','destination','time','geometry'])
    complete_routes_gdf = gpd.GeoDataFrame(complete_routes_df, crs = 'epsg:4326')
        
    return [complete_count_gdf, complete_routes_gdf]


# In[43]:


#results = generate_complete_edges_and_routes(OD_df[4784044133])


# In[44]:


OD_df_dict


# In[46]:


for key in OD_df_dict:
    print(key)


# In[ ]:


results = generate_complete_edges_and_routes(OD_df_dict[3935302581])


# In[ ]:


#import time

#start = time.time()

#results = {}

#count = 0
#for key in OD_df_dict:
    #while count > 1 and count < 6:
    #results[key] = generate_complete_edges_and_routes(OD_df_dict[key])
    #count += 1
    
#print(time.time() - start)


# In[ ]:


#for key in results:
# print edges
#print(results[key][0])
file_name = "./airport_output_edges_20200628_traffic_mean/weighted_sec_saved_edges_3935302581_traffic_mean.shp"
print('print file name:')
print(file_name)
results[0].to_file(driver = 'ESRI Shapefile', filename = file_name)


# In[ ]:


#results[1243386867][0]


# In[ ]:


# save as shapefile
#results[0].to_file(driver = 'ESRI Shapefile', filename = "./weighted_sec_saved_edges.shp")


# In[ ]:


# save as shapefile
#results[1].to_file(driver = 'ESRI Shapefile', filename = "./shortest_path_routes.shp")

