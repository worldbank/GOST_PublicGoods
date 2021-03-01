### other libraries

import pandas as pd

from shapely.ops import transform
from shapely.geometry import mapping, Polygon, box, shape
import matplotlib.pyplot as plt

from collections import namedtuple

from gbdxtools import Interface
from gbdxtools.task import env
from gbdxtools import CatalogImage

gbdx = Interface()

# input is a shapely polygon of the park of interest
def get_OSM_Amenities(parkid,park_shape):
    
    park_wkt=park_shape.wkt
    
#    Create a dataframe to add the counted amenities  
    columns=["Park_Id","Fac_Bench","Fac_Waste","Fac_Toilet","Fac_Water","Fac_Play","Fac_Hist","Fac_Retail","Fac_Fountain","Fac_Sports"]
    amenities_df = pd.DataFrame(columns=columns)

# query OSM vectors (results come back formatted as geojson)
    geojson = gbdx.vectors.query(park_wkt, query="ingest_source:OSM AND item_type:*", index="vector-osm-*",count=1E6)

#check if the amenities are within the park
    geom_list = []
    info_list = []
    for geojson in geojson:
        geom = shape(geojson['geometry'])
        info =geojson['properties']
        if park_shape.contains(geom)==True:
            info_list.append(info)
            
# Filter so there is only a list of item types left in the list
    types_list = []      
    for i in info_list:
        itemt= i[u'item_type']
        types_list.append(itemt)
        
# flatten the list so there is no list in list
    flat_list = [item for sublist in types_list for item in sublist]

# Count the number of amenities within each category
    Park_Id = parkid
    Fac_Bench = flat_list.count(u'Bench') + flat_list.count(u'Picnic Table')
    Fac_Waste = flat_list.count(u'Waste Basket')
    Fac_Toilet = flat_list.count(u'Toilets')
    Fac_Water = flat_list.count(u'Water Point') + flat_list.count(u'Drinking Water')
    Fac_Play = flat_list.count(u'Playground')
    Fac_Hist = flat_list.count(u'Historic (Memorial)') + flat_list.count(u'Historic (Monument)')
    Fac_Retail = flat_list.count(u'Retail') + flat_list.count(u'Cafe') + flat_list.count(u'Restaurant')
    Fac_Fount = flat_list.count(u'Fountain')
    Fac_Sports =  flat_list.count(u'Pitch')
    
# Append all counts to the dataframe
    amenities_df = amenities_df.append({"Park_Id":Park_Id,"Fac_Bench":Fac_Bench,"Fac_Waste":Fac_Waste,"Fac_Toilet":Fac_Toilet,"Fac_Water":Fac_Water,"Fac_Play":Fac_Play,"Fac_Hist":Fac_Hist,"Fac_Retail":Fac_Retail,"Fac_Fountain":Fac_Fount,"Fac_Sports":Fac_Sports},ignore_index=True)


    return amenities_df