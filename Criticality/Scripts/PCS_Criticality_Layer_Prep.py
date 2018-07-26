# -*- coding: utf-8 -*-
###################################################################################################
# Generate automatically layers of interest for use in criticality.
# Charles Fox October 2017
# Purpose: determine the criticality score for each linear feature in the network dataset
###################################################################################################
import os, sys, inspect, logging
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.geometry.base
import shapely.wkt
path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
path = os.path.split(path)[0]
dash = os.path.join(path,r'dashboard.xlsm')

ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)
critdf = pd.read_excel(dash, sheetname = "CRITICALITY", index_col = 0)
district = ctrl['Weight'].loc['DISTRICT']

logging.basicConfig(filename = os.path.join(path, 'PCS\Criticality\Runtime', "PCS_Criticality_Layer_Prep.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
logging.info("Starting Layer Prep Process")
print "Running: Criticality Layer Preparation Tool for %s. Do not interrupt" % district

#Settings
crs_in = {'init': 'epsg:4326'}   #WGS 84

#Folders
OD_IN = os.path.join(path, 'PCS\Criticality\Input', '%s' % district)
DATA_IN = os.path.join(path, 'PCS\Criticality\Data_Layers')
NETWORK_IN = os.path.join(OD_IN, 'Network.csv')
if not os.path.isdir(OD_IN):
        os.mkdir(OD_IN)
#Files
inNetwork = pd.read_csv(NETWORK_IN)
ginNetwork = gpd.GeoDataFrame(inNetwork,crs = crs_in, geometry = inNetwork['Line_Geometry'].map(shapely.wkt.loads))
inAdmin = os.path.join(DATA_IN,'Poverty_Communes_2009.shp')

for curFile in [dash, path, inAdmin, NETWORK_IN, DATA_IN,OD_IN,NETWORK_IN]:
    if not os.path.exists(curFile):
        logging.error("No input found: %s" % curFile)
        raise ValueError("No input found: %s" % curFile)

SingleNetworkObj = pd.DataFrame([str(ginNetwork.unary_union)], columns = ['geom'])
gSingleNetworkObj = gpd.GeoDataFrame(SingleNetworkObj,crs = crs_in, geometry = SingleNetworkObj['geom'].map(shapely.wkt.loads))
ginAdmin = gpd.read_file(inAdmin)
ginAdmin = ginAdmin.to_crs(crs_in)
ginAdmin = ginAdmin[['P_EName','D_EName','EN_name','geometry']]
ginAdmin['ID'] = ginAdmin.index
SelectedAdmins = gpd.sjoin(gSingleNetworkObj, ginAdmin, how="inner",op='intersects')
SelectedAdmins = SelectedAdmins.drop(['geom','geometry'], axis = 1)
ginAdmin = ginAdmin.loc[ginAdmin['ID'].isin(SelectedAdmins['ID']) == True]

def PrepareLayer(ginAdmin, layername, OD_IN):
    layer = os.path.join(DATA_IN, '%s.shp' % layername)
    glayer = gpd.read_file(layer)
    glayer = glayer.to_crs(crs_in)
    glayer['IDglayer'] = glayer.index
    SelectedFeatures = gpd.sjoin(ginAdmin, glayer, how="inner",op='intersects')
    glayer = glayer.loc[glayer['IDglayer'].isin(SelectedFeatures['IDglayer']) == True]
    glayer['Default'] = 1
    out = os.path.join(OD_IN,'%s.shp' % layername)
    glayer.to_file(out, driver = 'ESRI Shapefile')
    logging.info('Successfully generated local layer for %s' % layername)

for z in critdf.index:
    if (((critdf['DefaultLayer'][z]) != 0) & ((critdf['PrepDefault'][z]) == 'yes')):
        try:
            PrepareLayer(ginAdmin, critdf['DefaultLayer'][z], OD_IN)
        except:
            logging.info('no intersecting features for layer "%s" this road network' % critdf['DefaultLayer'][z])
            pass
