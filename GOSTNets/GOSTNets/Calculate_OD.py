import os, sys, logging, warnings, time

import osmnx
import networkx as nx
import pandas as pd
import geopandas as gpd
import numpy as np

from shapely.geometry import Point

import GOSTnet as gn

def calculateOD_gdf(G, origins, destinations, fail_value=-1, weight="time"):
    ''' Calculate Origin destination matrix from GeoDataframes
    
    Args:
        G (networkx graph): describes the road network. Often extracted using OSMNX
        origins (geopandas dataframe): source locations for calculating access
        destinations (geopandas dataframe): destination locations for calculating access
    Returns:
        numpy array: 2d OD matrix with columns as index of origins and rows as index of destinations
    '''
    #Get a list of originNodes and destinationNodes
    wgs84 = {'init':'epsg:4326'}
    if origins.crs != wgs84:
        origins = origins.to_crs(wgs84)
    if destinations.crs != wgs84:
        destinations = destinations.to_crs(wgs84)
    origins = gn.pandana_snap(G, origins)
    destinations = gn.pandana_snap(G, destinations)
    oNodes = origins['NN'].unique()
    dNodes = destinations['NN'].unique()
    od = gn.calculate_OD(G, oNodes, dNodes, fail_value)
    #TODO: the OD matrix is related to the unique list of origins and destinations
    #   instead of the origins and destination indices. How to fix that?
    # HERE'S HOW MUTHERFUCKER
    origins['OD_O'] = origins['NN'].apply(lambda x: np.where(oNodes==x)[0][0])
    destinations['OD_D'] = destinations['NN'].apply(lambda x: np.where(dNodes==x)[0][0])
    return(od[origins['OD_O'].values,:][:,destinations['OD_D'].values])    
    
def calculateOD_csv(G, originCSV, destinationCSV='', 
        oLat="Lat", oLon="Lon", dLat="Lat", dLon="Lon",
        crs={'init':'epsg:4326'},
        fail_value=-1, weight='time'):
    ''' Calculate OD matrix from csv files of points
    Args:
        G (networkx graph): describes the road network. Often extracted using OSMNX
        origins (string): path to csv file with locations for calculating access
        destinations (string): path to csv with destination locations for calculating access
        oLat/oLon/dLat/dLon (string, optional): name of latitude and longitude in origin/destination files.
            defaults to "Lat" and "Lon for both origins and destinations
        crs (dictionary, optional): crs of input origins and destinations, defaults to {'init':'epsg:4326'}
        fail-value (integer, optional): value to put in OD matrix if no route found, defaults to -1
        weight (string, optional): variable in G used to define edge impedance, defaults to 'time'
    Returns:
        numpy array: 2d OD matrix with columns as index of origins and rows as index of destinations
    '''
    originPts = pd.read_csv(originCSV)
    pts = [Point(x) for x in zip(originPts[oLon],originPts[oLat])]
    originGDF = gpd.GeoDataFrame(originPts, geometry=pts, crs=crs)
    if destinationCSV == '':
        destinationGDF = originGDF.copy()
    else:
        originPts = pd.read_csv(destinationCSV)
        pts = [Point(x) for x in zip(originPts[dLon],originPts[dLat])]
        destinationGDF = gpd.GeoDataFrame(originPts, geometry=pts, crs=crs)
    OD = calculateOD_gdf(G, originGDF, destinationGDF, fail_value, weight)
    return(OD)
   
def calculate_gravity(od, oWeight=[], dWeight=[], decayVals=[0.01,
                                                        0.005,
                                                        0.001,
                                                        0.0007701635,   # Market access halves every 15 mins
                                                        0.0003850818,   # Market access halves every 30 mins
                                                        0.0001925409,   # Market access halves every 60 mins
                                                        0.0000962704,   # Market access halves every 120 mins
                                                        0.0000385082,   # Market access halves every 300 mins
                                                        0.00001]):
    ''' Calculate the gravity weight between origins and destinations.
    
    Args:
        od (ndarray): matrix of travel time    
        oWeight/dWeight (array, optional) - array of weights for calculating weight; reverts to 1 for if not defined   
        decayVals (array, optional): decayVals to calculate for gravity. Each value will be returned as a column of results
    Returns:
        geopandas: columns of decayvals
    '''
    if len(oWeight) != od.shape[0]:
        oWeight = [1] * od.shape[0]
    if len(dWeight) != od.shape[1]:
        dWeight = [1] * od.shape[1]
    allRes = []
    for dist_decay in decayVals:
        outOD = od * 0
        decayFunction = lambda x: np.exp(-1 * dist_decay * x)
        for row in range(0, od.shape[0]):
            curRow = od[row,:]
            decayedRow = decayFunction(curRow)
            weightedRow = decayedRow * oWeight[row] * dWeight
            outOD[row,:] = weightedRow
        summedVals = np.sum(outOD, axis=1)
        allRes.append(summedVals)
    res = pd.DataFrame(allRes).transpose()
    res.columns = columns=['d_%s' % d for d in decayVals]
    return(res)
    
    