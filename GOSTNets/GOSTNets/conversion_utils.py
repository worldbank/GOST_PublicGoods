import os, sys
import rasterio

import pandas as pd
import geopandas as gpd
import numpy as np

from scipy import interpolate
from rasterio import features
from rasterio.mask import mask
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling

def rasterize_od_results(inD, outFile, field):
    ''' Convert gridded point data frame to raster of commensurate size and resolution
    
    INPUT
    inD [ geopandas data frame ] - OD matrix as point data frame
    outFile [ string ] - path to save output raster
    field [ string ] - field to rasterize
    
    RETURNS
    None
    '''
    #create grid from input shapefile
    # get xs, ys, and values from origin points
    xs = np.array(inD.geometry.apply(lambda p: p.x))
    ys = np.array(inD.geometry.apply(lambda p: p.y))
    vals = np.array(inD[field])
    # creates a full grid for the entire bounding box (all pairs of xs and ys)
    unique_xs = np.unique(xs)
    unique_ys = np.unique(ys)
    xx, yy = np.meshgrid(unique_xs, unique_ys)
    # this creates a new set of values to fill the grid
    grid_array = interpolate.griddata((xs,ys), vals, (xx, yy))
    x_pixels = grid_array.shape[1]
    y_pixels = grid_array.shape[0]
    # get the right transformation for raster file
    xRes = (xx.max() - xx.min()) / len(unique_xs)
    yRes = (yy.max() - yy.min()) / len(unique_ys)
    # get the right transformation for raster file
    trans = rasterio.transform.from_bounds(xx.min() - (xRes/2), yy.min() - (yRes/2), 
                                           xx.max() - (xRes/2), yy.max() - (yRes/2),
                                           x_pixels - 1, y_pixels - 1)
    new_dataset = rasterio.open(
        outFile, 'w', driver = 'GTiff',
        height = y_pixels, width = x_pixels,
        count=1, dtype=str(grid_array.dtype),
        crs=inD.crs,
        transform=trans
    )
    
    shapes = ((row.geometry,row[field]) for idx, row in inD.iterrows())
    burned = features.rasterize(shapes=shapes, fill=0, out_shape=grid_array.shape, transform=new_dataset.transform)
    burned = burned.astype(grid_array.dtype)    
    new_dataset.write_band(1, burned)
    
    new_dataset.close()