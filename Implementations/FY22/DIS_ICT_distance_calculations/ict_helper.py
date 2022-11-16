import sys, os, importlib, json, multiprocessing, time
import rasterio, pycountry
import reverse_geocode

import geopandas as gpd
import pandas as pd

from urllib.request import urlopen
from shapely.geometry import Point
from shapely.ops import nearest_points
from shapely import wkt

# Import GOST libraries; sys.path.append will be unnecessary if libraries are already installed
sys.path.append("../../../../gostrocks/src")
sys.path.append("../../../../GOST_Urban/src")

import GOSTRocks.rasterMisc as rMisc
from GOSTRocks.misc import tPrint
import GOST_Urban.UrbanRaster as urban

def get_nearest_date(x, selCables):
    selCables = selCables.loc[selCables['RFS'] <= x['d2_l1_year_perf_indicators']]
    if selCables.shape[0] > 0:
        xx = selCables.loc[selCables.geometry == nearest_points(x['geometry'], selCables.unary_union)[1]]
        return(xx.sort_values(['RFS'])['RFS'].iloc[0])
    else:
        return(-1)

def get_nearest(x, selCables):
    selCables = selCables.loc[selCables['RFS'] <= x['d2_l1_year_perf_indicators']] 
    if selCables.shape[0] == 0:
        return(-1)
    return(x['geometry'].distance(selCables.unary_union))

def calculate_country(iso3, curD, curB, curN, out_file, selCol, selIXP, inCables, inCell, epsg='epsg:6933', debug=False):
    ''' calculaet ICT distances per country
    
    Args:
        curD: geopandas data frame of WBES survey locations
        curB: geopandas data frame of country bounds        
        curN: geopandas data frame of neighbouring countries boundaries
        outFile: string of the path for the output file; is read in if it doesn't exist
        selCol: geopandas data frame of colocation centers 
        selIXP: geopandas data frame of IXPs 
        inCables: geopandas data frame cable landing spots 
    '''
    start_time = time.time()
    tPrint("Starting %s" % iso3)
    cell_coverage_folder = '/home/public/Data/GLOBAL/INFRA/GSMA/2019/MCE/Data_MCE/Global'
    cell_files = ['MCE_Global2G_2020.tif', 'MCE_Global3G_2020.tif', 'MCE_Global4G_2020.tif']
    gsma2g_R = rasterio.open(os.path.join(cell_coverage_folder, cell_files[0]))
    gsma3g_R = rasterio.open(os.path.join(cell_coverage_folder, cell_files[1]))
    gsma4g_R = rasterio.open(os.path.join(cell_coverage_folder, cell_files[2]))
    if False: #os.path.exists(out_file):
        curD = pd.read_csv(out_file, index_col=0)
        curD = pd.merge(distD, curD.loc[:,['idstd','d2_l1_year_perf_indicators']], on='idstd')
        curD_geom = curD['geometry'].apply(wkt.loads)
        distD = gpd.GeoDataFrame(curD, geometry=curD_geom, crs=epsg)
        # Remove columns that need to be re-calculated
        distD = distD.loc[:,[not "ngh" in x for x in distD.columns]]
        distD = distD.loc[:,[not "gsma" in x for x in distD.columns]]
        distD = distD.loc[:,[not "cables_dist" in x for x in distD.columns]]       
    else:
        distD = curD.to_crs(epsg)
        
    total_bound = curB.unary_union
    if curB.shape[0] > 0:
        if not 'col_dist' in distD.columns:
            if selCol.shape[0] > 0:
                selCol = selCol.to_crs(epsg)
                distD['col_dist'] = distD.distance(selCol.unary_union)
            else:
                distD['col_dist'] = -1

        if not "ixp_dist" in distD.columns:
            if selIXP.shape[0] > 0:
                selIXP = selIXP.to_crs(epsg)
                distD['ixp_dist'] = distD.distance(selIXP.unary_union)
            else:
                distD['ixp_dist'] = -1

        if not 'firstCable' in distD.columns:
            selCables = inCables.loc[inCables['ISO3'] == iso3]            
            if selCables.shape[0] > 0:
                selCables = selCables.to_crs(epsg)
                # Calculate distance and date of first cable landing point
                first_date = selCables['RFS'].sort_values().iloc[0]
                first_points = selCables.loc[selCables['RFS'] == first_date]
                distD['firstCable'] = first_date
                distD['firstCable_dist'] = distD.distance(first_points.unary_union)
                # Calculate distance and date of closest cable landing point
                distD['closestCable'] = distD.apply(lambda x: get_nearest_date(x, selCables), axis=1)                
                distD['closestCable_dist'] = distD.apply(lambda x: get_nearest(x, selCables), axis=1)
            else:
                distD['firstCable'] = ''
                distD['firstCable_dist'] = -1
                # Calculate distance and date of closest cable landing point
                distD['closestCable'] = ''
                distD['closestCable_dist'] = -1

        # Calculate distance to nearest neighbouring country
        if not "ngh1_dist" in distD.columns:
            cnt = 1
            for idx, row in curN.iterrows():                 
                distD['ngh%s' % cnt] = row['ISO3']
                distD['ngh%s_dist' % cnt] = distD.distance(row['geometry'])    
                #Calculate distance to neighbouring submarine cables
                selCables = inCables.loc[inCables['ISO3'] == row['ISO3']]
                if debug:
                    tPrint(f'Neighbour cable distance {row["ISO3"]}: {selCables.shape[0]}')
                if selCables.shape[0] > 0:
                    #distD['ngh%s_cbl_dist' % cnt] = distD.distance(selCables.unary_union)
                    distD['ngh%s_cbl_dist' % cnt] = distD.apply(lambda x: get_nearest(x, selCables), axis=1)
                    distD['ngh%s_cbl' % cnt] = distD.apply(lambda x: get_nearest_date(x, selCables), axis=1)
                else:
                    distD['ngh%s_cbl_dist' % cnt] = -1
                    distD['ngh%s_cbl' % cnt] = -1
                cnt = cnt +1            
        
        if not False: # 'cell_dist' in distD.columns:
            cell_sindex = inCell.sindex
            potential_matches = inCell.loc[list(cell_sindex.intersection(total_bound.bounds))]
            selCell = potential_matches.loc[potential_matches.intersects(total_bound)]
            selCell = selCell.to_crs(epsg)
            distD['cell_dist'] = distD.distance(selCell.unary_union)

        if not "gsma2g" in distD.columns:
            coordsD = distD.to_crs(gsma2g_R.crs)
            coords = [[x.x,x.y] for x in coordsD['geometry']]
            distD['gsma2g'] = [x[0] for x in list(gsma2g_R.sample(coords))]
            distD['gsma3g'] = [x[0] for x in list(gsma3g_R.sample(coords))]
            distD['gsma4g'] = [x[0] for x in list(gsma4g_R.sample(coords))]

        pd.DataFrame(distD).to_csv(out_file)
        return(distD)
    end_time = time.time()
    tPrint(f"Completed {iso3}: {round((end_time-start_time)/60)}")