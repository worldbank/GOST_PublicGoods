import sys, os, importlib, json, boto3
import rasterio, geojson, folium, h3

import pandas as pd
import geopandas as gpd

from shapely.geometry import Point, mapping

sys.path.insert(0, "/home/wb411133/Code/gostrocks/src")
import GOSTRocks.rasterMisc as rMisc
import GOSTRocks.ntlMisc as ntl
from GOSTRocks.misc import tPrint

#class wb_survey_handler():


class survey_extract():
    ''' Extract geospatial variables for survey locations 
    '''
    def __init__(self, inD, lat_col, lon_col, id_column, crs=4326):
        ''' Extract geospatial information for survey locations
        
        Args:
            inD: pandas dataframe of survey data
            lat_column: string name of column in survey_data containing latitude
            lon_column: string name of column in survey_data containing longitude
        '''
        self.surveyD = inD
        self.id_column = id_column
        # create GeoDataFrame of surveyD
        geoms = [Point(x) for x in zip(inD[lon_col], inD[lat_col])]
        inG = gpd.GeoDataFrame(inD, geometry=geoms, crs=crs)
        self.inG = inG        
        # Create coordinates list for use in rasterio.sample
        self.coords = [(x,y) for x,y in zip(inG.geometry.x, inG.geometry.y)]
        
    def calc_zonal(self, raster_file_defs, verbose=False):
        ''' extract point location values for raster files
        
        Args:
            raster_file_defs: list of lists; raster file path, out_column_name. Out_column_name 
                will be appended to input survey data to create results
        Returns:
            pandas data frame with survey id column and resulting out_column_names
        '''
        outD = self.surveyD.copy()
        out_cols = [self.id_column]
        for file_def in raster_file_defs:
            inR = rasterio.open(file_def[0])
            outD[file_def[1]] = [x[0] for x in inR.sample(self.coords)]
            out_cols.append(file_def[1])
            if verbose: 
                tPrint(file_def[1])
            
            
        outD = outD.loc[:,out_cols]
        return(outD)
        
    