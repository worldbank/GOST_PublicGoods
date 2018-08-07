################################################################################
# Process routing via OSRM
# Benjamin Stewart, August 2018
# Purpose: Calculate a number of OSRM tasks 
################################################################################

import json, sys, os, time, argparse
import shapely, geojson

import pandas as pd
import geopandas as gpd
import numpy as np
import urllib2 as url  

from shapely.geometry import shape

class osrm_service:
    def __init__(self, **kwargs):
        if 'osrmHeader' in kwargs.keys():
            self.osrmHeader = kwargs['osrmHeader']
        else:
            self.osrmHeader = 'http://router.project-osrm.org'
            
    def prepMatrix(self, matrixDF, oDF, dDF):
        ''' Create a DataFrame for use in extractGeoms from a OD matrix
        --- variables ---
        matrixDF [pandas Dataframe] - input OD matrix
        oDF [pandas Dataframe] - details of source features from dataframe
        dDF [pandas Dataframe] - details of dest   features from dataframe
        --- returns ---
        [pandas DF]
        '''
        # remove 0 distances from matrix
        matrixDF = matrixDF[matrixDF['DIST'] != 0]
        # df.loc[df.groupby("item")["diff"].idxmin()]
        minDF = matrixDF.loc[matrixDF.groupby('O_UID')['DIST'].idxmin()]

        minDF = minDF.merge(oDF, left_on="O_UID", right_on="ID")
        minDF = minDF.rename(columns = {'Lat':'sLat','Lon':'sLon'})

        minDF = minDF.merge(dDF, left_on="D_UID", right_on="ID")
        minDF = minDF.rename(columns = {'Lat':'dLat','Lon':'dLon'})
        return minDF           
    
    def extractGeoms(self, pairsDF, sLat="sLat", sLon="sLon", dLat="dLat", dLon="dLon", verbose=False):
        ''' Runs the routing command on all the rows in the inputDF        
        --- variables ---
        pairsDF [pandas DataFrame] - dataframe containing two sets of latitude and 
            longitude as defined by the optional commands
        [optional] sLat ... [string] - field name containing the lats and longs
        --- returns ---
        [geopandas dataframe] - contains two(ish?) columns - the index of the input pairs DF
        '''
        allRes = []
        for idx, row in pairsDF.iterrows():
            if verbose:
                print idx
            # 94.03333,19.78333;108.212426,16.051264?overview=simplified&geometries=geojson
            location = "%s,%s;%s,%s" % (row[sLon],row[sLat],row[dLon],row[dLat])
            baseUrl = "%s/route/v1/driving/%s?overview=simplified&geometries=geojson" % (self.osrmHeader, location)
            try:
                r = url.urlopen(baseUrl)
                curRes = json.loads(r.read().decode('utf-8'))
                curCoords = curRes['routes'][0]['geometry']['coordinates']
                curGeom = geojson.LineString(curCoords)
                allRes.append([idx, shape(curGeom)])
            except:
                print("There was an error with row %s" % idx)
        final = pd.DataFrame(allRes, columns=['ID','geometry'])
        final = gpd.GeoDataFrame(final, geometry='geometry')
        return final

if __name__ == "__main__":
    xx = osrm_service()
    
    inFolder = r"C:\Users\WB411133\OneDrive - WBG\AAA_BPS\Code\Code\Github\GOST_PublicGoods\Market_Access"
    inOd = pd.read_csv(os.path.join(inFolder, "Os.csv"))
    inDd = pd.read_csv(os.path.join(inFolder, "Ds.csv"))
    inMatrix = pd.read_csv(os.path.join(inFolder, "OD_Matrix.csv"))
    
    inRes = xx.prepMatrix(inMatrix, inOd, inDd)        
    curRoute = xx.extractGeoms(inRes)
    print curRoute
    curRoute.to_csv(os.path.join(inFolder, "nearest_paths.csv"))