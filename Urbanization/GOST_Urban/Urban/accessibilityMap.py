################################################################################
# Calculate Accessibility using Global Data
# Benjamin P. Stewart, June 2017
# Purpose: Calculate accessibility from Fay, Deichmann, etc. based on network distance
#   and population of market access points.
################################################################################

import os, sys
import shapefile
import osrm
import pandas

class accessMap(object):
    '''
    An accessibility map describes a population weighted measure of access to markets.
    This is based on a combination of origins and destinations, and the network that 
    connects them. 
    
    The default destinations are identical to the origins, in order to measure access 
    between all markets, but can be modified to measure other kinds of access. The default
    network access methodology uses OSRM and the OSM network (so use with caution!)
    '''
    
    def __init__(self, origins, destinations=-1):
        '''
        origins - the input pandas data from containing the point-based market locations.
            The data frame should be organized as three columns [ID, Long, Lat, Weight]
            --> Weight - should be population, but could also be other measures of gravity 
                              for the destinations shapefile
        destinations - if set to -1, destinations match origins, otherwise should match formatting
            of the origins dataset
        '''
        self.origins = origins
        if destinations == -1:
            self.destinations = origins
        else:
            self.destinations = destinations
            
    def calcODMatrix(self, longName="Long", latName="Lat", idName="id"):
        '''
        Generates the OSRM table based on the origins and destinations. The x and y locations
        need to be extracted from the shapefiles
        longName - column name of longitude
        latName - column name of latitude
        idName - column name of id
        '''
        list_coord = [tuple(x) for x in self.origins[[longName, latName]].values]
        list_id = self.origins[[idName]].values
        
        print list_coord
        print list_id
        time_matrix, snapped_coords = osrm.table(list_coord, ids_origin=list_id, output="dataframe")
        
        self.od_matrix = time_matrix
        self.snappedOrigins = snapped_coords
        return self.odMatrix
    
    def accessibilityTable(self, newWeights=-1):
        '''
        based on the odMatrix, calculate the accessibility matrix
        weights - if not defined, the weights come from the origins dataframe
        '''
        if not self.od_matrix:
            odMatrix = calcODMatrix(self)
        else:
            odMatrix = self.od_matrix
        if newWeights != -1:
            weights = newWeights
        else:
            weights = self.origins[[self.origins.columns.values[3]]]
        
        