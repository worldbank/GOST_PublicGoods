import sys, os, importlib
import fiona

import pandas as pd
import geopandas as gpd
import numpy as np

import GOSTRocks.rasterMisc as rMisc
import GOSTRocks.misc as misc

from shapely.geometry import Point, box, Polygon
from shapely.wkt import loads
from math import ceil


def extract_da_buildings(inAOI, in_da, out_crs):
    ''' Extract buildings from digitize africa buildings dataset for specific AOI
    
    INPUT
        inAOI [geopandas dataframe] - AOI of interest
        in_da [ string ] - path to shapefile containing DA buildings
        out_crs [ string ] - string representation of output CRS as 'epsg:3857'
        
    RETURNS
    '''
    all_da = gpd.read_file(in_da)
    allda_idx = all_da.sindex
    if inAOI.crs != all_da.crs:
        inAOI = inAOI.to_crs(all_da.crs)
    potential_da = all_da.loc[list(allda_idx.intersection(inAOI.total_bounds))]
    sel_da = potential_da.loc[potential_da['geometry'].apply(lambda x: x.intersects (inAOI.unary_union))]
    sel_da = sel_da.to_crs(out_crs)
    return(sel_da)
    
def generate_grid(inAOI, res, out_crs):
    ''' generate a vector grid of resolution RES intersecting features in inAOI
    
        INPUT
            inAOI[ geopandas dataframe ] - AOI
            res [ number ] - resolution of output grid in unit of out_crs
            out_crs [ string ] - string representation of output CRS as 'epsg:3857'
        RETURNS
            [ geopandas dataframe ] - vector grid of res RES intersecting inAOI        
    '''
    if str(inAOI.crs) != out_crs:
        inAOI = inAOI.to_crs(out_crs)
    xmin, ymin, xmax, ymax = inAOI.total_bounds
    gridWidth, gridHeight = [res, res]
    
    cols = list(np.arange(xmin, xmax + gridWidth, gridWidth))
    rows = list(np.arange(ymin, ymax + gridHeight, gridHeight))
    all_res = []
    
    allP = inAOI.unary_union
    
    for x in cols:
        for y in rows:
            poly = Polygon([(x,y), (x+gridWidth, y), (x+gridWidth, y+gridHeight), (x, y+gridHeight)])
            if poly.intersects(allP):
                all_res.append([x,y,poly])        
    grid = gpd.GeoDataFrame(pd.DataFrame(all_res, columns=['rowIdx', 'colIdx', 'geometry']), geometry='geometry', crs=out_crs)
    return(grid)
    
def summarize_in_grid(grid, inDA, inB, inP = None):
    ''' summarize the locations of buildings and parcels inside a regularly spaced grid
    
    INPUT
        grid [ geopandas dataframe ] - generated from helper function generate_grid
        inDA [ geopandas dataframe ] - buildings from Digitize Africa
        inB [ geopandas dataframe ] - buildings from collected source, can be points or polygons
        inP [ geopandas dataframe ] - collected parcels
    '''
    # Check projections
    if inB.crs.to_epsg() != grid.crs.to_epsg():
        raise(ValueError("inB CRS must match grid"))
    if inDA.crs.to_epsg() != grid.crs.to_epsg():
        raise(ValueError("inDA CRS must match grid"))
    
    try:
        grid.reset_index(inplace=True)
        inDA.reset_index(inplace=True)
        inB.reset_index(inplace=True)
    except:
        pass
        
    grid['per_b'] = 0.
    grid['per_da']= 0.
    if not inP is None:
        grid['per_p'] = 0.
        inP.reset_index(inplace=True)
        p_idx = inP.sindex
    
    da_idx = inDA.sindex
    b_idx = inB.sindex
    
    b_type = inB.geom_type.iloc[0]
    if b_type == "Point":
        inDA['geometry'] = inDA['geometry'].apply(lambda x: x.centroid)
        
    for idx, row in grid.iterrows():
        # identify intersecting da buildings 
        potential_da = inDA.loc[list(da_idx.intersection(row['geometry'].bounds))]    
        i_da = potential_da.loc[potential_da.intersects(row['geometry'])]
        c_da = potential_da.loc[potential_da['geometry'].apply(lambda x: row['geometry'].contains(x))]
                    
        # identify intersecting buildings 
        potential_buildings = inB.loc[list(b_idx.intersection(row['geometry'].bounds))]    
        i_bld = potential_buildings.loc[potential_buildings.intersects(row['geometry'])]
        c_bld = potential_buildings.loc[potential_buildings['geometry'].apply(lambda x: row['geometry'].contains(x))]
        
        if not inP is None:
            # identify intersecting parcels
            potential_parcels = inP.loc[list(p_idx.intersection(row['geometry'].bounds))]    
            i_par = potential_parcels.loc[potential_parcels.intersects(row['geometry'])]
            c_par = potential_parcels.loc[potential_parcels['geometry'].apply(lambda x: row['geometry'].contains(x))]
            
        # calulate percent parcel and percent built
        try:
            if b_type == "Point":
                per_da = i_da.shape[0]
            else:
                per_da = row['geometry'].intersection(i_da.unary_union).area/row['geometry'].area
        except:
            per_da = 0
        
        try:
            if b_type == "Point":
                per_b = i_bld.shape[0]
            else:
                per_b = row['geometry'].intersection(i_bld.unary_union).area/row['geometry'].area
        except:
            per_b = 0
        try:
            per_parcel   = row['geometry'].intersection(i_par.unary_union).area/row['geometry'].area
        except:
            per_parcel = 0
            
        grid.loc[idx, 'per_b'] = per_b
        if not inP is None:
            grid.loc[idx, 'per_p'] = per_parcel
        grid.loc[idx, 'per_da']= per_da
    return(grid)
        
def compare_buildings_parcels(inP, inB, inDA, main_size = 50):
    inB['area'] = inB['geometry'].apply(lambda x: x.area)
    inDA['area'] = inDA['geometry'].apply(lambda x: x.area)
    # Summarize buildings in parcels
    inP['BLDG_I'] = 0 # buildings intersecting parcel
    inP['BLDG_C'] = 0 # buildings contained in parcel
    inP['BLDG_M'] = 0 # buildings in parcel > 50m2

    inP['BLDG_I_DA'] = 0 # same as above, but from digitize africa dataset
    inP['BLDG_C_DA'] = 0
    inP['BLDG_M_DA'] = 0
    
    b_idx = inB.sindex
    g_idx = inDA.sindex
    for idx, row in inP.iterrows():    
        # Summarize collected buildings
        potential_buildings = inB.loc[list(b_idx.intersection(row['geometry'].bounds))]
        m_bld = potential_buildings.loc[~potential_buildings.intersects(row['geometry'])]
        i_bld = potential_buildings.loc[potential_buildings.intersects(row['geometry'])]
        c_bld = potential_buildings.loc[potential_buildings['geometry'].apply(lambda x: row['geometry'].contains(x))]    
        inP.loc[idx, 'BLDG_I'] = i_bld.shape[0]
        inP.loc[idx, 'BLDG_C'] = c_bld.shape[0]
        m_bld = i_bld.loc[i_bld['area'] > 50]
        inP.loc[idx, 'BLDG_M'] = m_bld.shape[0]
        
        # Summarize buildings in Google
        potential_buildings = inG.loc[list(g_idx.intersection(row['geometry'].bounds))]
        m_bld_g = potential_buildings.loc[~potential_buildings.intersects(row['geometry'])]
        i_bld_g = potential_buildings.loc[potential_buildings.intersects(row['geometry'])]
        c_bld_g = potential_buildings.loc[potential_buildings['geometry'].apply(lambda x: row['geometry'].contains(x))]    
        inP.loc[idx, 'BLDG_I_DA'] = i_bld_g.shape[0]
        inP.loc[idx, 'BLDG_C_DA'] = c_bld_g.shape[0]
        m_bld_g = i_bld_g.loc[i_bld_g['area'] > main_size]
        inP.loc[idx, 'BLDG_M_DA'] = m_bld_g.shape[0]
    return(inP)

def attribute_buildings_parcels(inB, inP):
    # Attribute buildings with parcel stats
    p_idx = inP.sindex
    inB['p_i'] = 0
    inB['p_c'] = 0
    max_idx = inB.index[-1]
    for idx, row in inB.iterrows():
        print(f'{idx} of {max_idx}')
        potential_parcels = inP.loc[list(p_idx.intersection(row['geometry'].bounds))]
        i_par = potential_parcels.loc[potential_parcels.intersects(row['geometry'])]
        c_par = potential_parcels.loc[potential_parcels['geometry'].apply(lambda x: x.contains(row['geometry']))]
        
        inB.loc[idx, 'p_i'] = i_par.shape[0]
        inB.loc[idx, 'p_c'] = c_par.shape[0]
    return(inB)