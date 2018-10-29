import sys, os

import osmnx as ox
import geopandas as gpd

sys.path.insert(0, r"C:\Users\WB411133\OneDrive - WBG\AAA_BPS\Code\Code\Github\GOST_PublicGoods")
from Market_Access import OSMNX_POIs

inD = gpd.read_file(r"Q:\WORKINGPROJECTS\CityScan\Data\CityExtents\Bamako_AOI.shp")
inD = inD.to_crs({'init': 'epsg:4326'})

current = OSMNX_POIs.AmenityObject('Education', inD.unary_union, ['school','university','secondary school', 'kindergarten', 'college'], 
                                    "C:\Temp")
            
temp = current.GenerateOSMPOIs()

relations = [{'type': 'relation', 'members': [{'type': 'way', 'role': 'outer', 'ref': 96565478}, {'type': 'way', 'role': 'inner', 'ref': 206575417}], 'tags': {'type': 'multipolygon', 'amenity': 'school', 'name': 'ex charl de gaul'}, 'id': 2776558}, {'type': 'relation', 'members': [{'type': 'way', 'role': 'outer', 'ref': 313403837}, {'type': 'way', 'role': 'outer', 'ref': 313403834}], 'tags': {'type': 'multipolygon', 'amenity': 'university', 'name': 'Université de Bamako'}, 'id': 4214679}, {'type': 'relation', 'members': [{'type': 'node', 'role': '', 'ref': 3916380819}], 'tags': {'amenity': 'kindergarten', 'source': 'survey', 'name': 'Creche Maternelle NFALY Sacko'}, 'id': 5808364}, {'type': 'relation', 'members': [{'type': 'way', 'role': 'outer', 'ref': 295305317}, {'type': 'way', 'role': 'inner', 'ref': 592707865}], 'tags': {'type': 'multipolygon', 'amenity': 'school', 'building': 'school', 'name': 'Complexe Scolaire Platon'}, 'id': 8344131}]