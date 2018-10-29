import os, sys, logging, json

import rasterio, geojson

import geopandas as gpd
import pandas as pd
import numpy as np

from rasterio.features import shapes
from rasterio.mask import mask
from shapely.geometry import shape

def getFeatures(gdf, idx=0):
    #Function to parse features from GeoDataFrame in such a manner that rasterio wants them
    return [json.loads(gdf.to_json())['features'][idx]['geometry']]

#Get all values built now
def processGhsl(curR, minVal, date, outT, ghslType="FULL"):
    ''' Convert GHSL to vector
    INPUT
    curR [rasterio] - GHSL raster object clipped to study area
    minVal [int] - value to extract from the dataset 3 (2014), 4(2000), 5(1990), 6(1975)
    date [string]- date to include in returned geodataframe
    outT [Affine Transform] - affine transformation for curR
    [optional] ghslType [string] - type of data conversion: FULL means all built area for that date, 
                                    PART means just that year . You can also pass an array of specific
                                    values to check
    RETURNS
    [geopandas data frame]
    
    DEBUG
    curR = out_img
    minVal = 3
    date = "d014"
    outT = out_transform    
    '''
    if ghslType == "FULL":
        bNow = curR >= minVal
    elif ghslType == "PART":
        bNow = curR == minVal
    elif type(ghslType) == type([]):
        for v in ghslType:
            try:
                bNow = bNow + (curR == v)
            except:
                bNow = curR == v
    else:
        raise ValueError("What the fuck did you pass in for ghslType!?")
    bNowShape = shapes(bNow.astype(np.int16), transform=outT)
    allGeoms = [shape(geojson.loads(json.dumps(xx[0]))) for xx in bNowShape if xx[1] != 0]
    curDF = gpd.GeoDataFrame(pd.DataFrame({"Date":[date] * len(allGeoms)}), geometry=allGeoms)    
    return(curDF)

def calculateLEI(gNew, gOld, dist=300):
    ''' Calculate Landscape Expansion index on input rows
    REFERENCE: Xiaoping Liu, Xia Li, Yimin Chen, Zhangzhi Tan, Shaoying Li, Bin Ai
    A new landscape index for quantifying urban expansion using multi-temporal remotely sensed data
    Landscape Ecology (2010) 25:671-682
    
    INPUT
    gNew [geopandas dataframe]
    gOld [geopandas dataframe]

    Returns [List of double] - LEI for each features in gNew
    
    DEBUG:
    gNew = new2014
    gOld = existing2000
    idx = 3890
    row = gNew.iloc[idx]
    cGeom = row.geometry
    cGeomBuffer = cGeom.buffer(dist).difference(cGeom)    
    oldShape = gOld.unary_union
    
    '''
    LEIValues = []
    for idx, row in gNew.iterrows():
        if idx % 50 == 0:
            print ("Processing %s of %s: %s" % (idx, gNew.shape[0], idx/gNew.shape[0]))
        cGeom = row.geometry
        #Buffer and clip the current geometry
        cGeomBuffer = cGeom.buffer(dist).difference(cGeom)
        #Select intersecting features from gOld
        try:
            curGold = gOld[gOld.intersects(cGeomBuffer)].unary_union
            if not curGold.is_valid:
                curGold = curGold.buffer(0)
            cGeomIntersection = cGeomBuffer.intersection(curGold)
            LEI = (cGeomIntersection.area)/(cGeomBuffer.area)
        except:
            LEI = 0
        LEIValues.append(LEI*100)
    return LEIValues

if __name__ == "__main__":
    calculateLEI = False
    inGHSL = "Q:/GLOBAL/POP&DEMO/GHS/BETA/FULL/MT/MT.vrt"
    inShp = r"Q:\WORKINGPROJECTS\CityScan\Data\CityExtents\Conotou_AOI.shp"
    inShp = r"C:\Temp\LEI_Cities\Cities footprint.shp"
    temp_outFolder = "C:/Temp/LEI_Cities"

    inR = rasterio.open(inGHSL)
    meta = inR.meta.copy()
    inD = gpd.read_file(inShp)
    inD = inD.to_crs(inR.crs)

    #Loop through input dataset
    for idx, row in inD.iterrows():
        print(idx)
        if idx >= 0:            
            outFolder = os.path.join(temp_outFolder, "%s" % idx)
            try:
                #os.mkdir(outFolder)
                coords = getFeatures(inD, idx)    
                out_img, out_transform = mask(inR, shapes=coords, crop=True)        
                #Compare urban expansions to existing built area for three comparisons
                #   new2014     vs existing2000
                #   new2000     vs existing1990
                #   new20002014 vs existing1990
                    
                #Create all baseline data
                meta.update({"height":out_img.shape[1], "width":out_img.shape[2],
                                         "transform":out_transform, "driver":"GTiff"})
                with rasterio.open(os.path.join(outFolder, "GHSL.tif"), 'w', **meta) as outFile:
                    outFile.write(out_img)
                    
                new2014 = processGhsl(out_img, 3, "2014", out_transform, "PART")
                new2014.crs = inR.crs
                new2014.to_file(os.path.join(outFolder, "new2014.shp"))
                
                existing2000 = processGhsl(out_img, 4, "2014", out_transform, [4,5,6,1])
                existing2000.crs = inR.crs
                existing2000.to_file(os.path.join(outFolder, "existing2000.shp"))
                
                new2000 = processGhsl(out_img, 4, "2014", out_transform, "PART")
                new2000.crs = inR.crs
                new2000.to_file(os.path.join(outFolder, "new2000.shp"))

                existing1990 = processGhsl(out_img, 4, "2014", out_transform, [4,5,6,1])
                existing1990.crs = inR.crs
                existing1990.to_file(os.path.join(outFolder, "existing1990.shp"))
                
                new0014 = processGhsl(out_img, 4, "2014", out_transform, [3,4])
                new0014.crs = inR.crs
                new0014.to_file(os.path.join(outFolder, "new0014.shp"))
                
                if calculateLEI:
                    lei_2014_2000     = calculateLEI(new2014, existing2000, dist=300)
                    #Attache LEI results to dataframe
                    new2014['LEI'] = pd.DataFrame(lei_2014_2000)
                    new2014.to_file(leiOutFile)
                    
                    lei_2000_1990     = calculateLEI(new2000, existing1990, dist=300)
                    #Attache LEI results to dataframe
                    new2000['LEI'] = pd.DataFrame(lei_2000_1990)
                    new2000.to_file(os.path.join(outFolder, "LEI_2000_1990.shp"))
                    
                    lei_20142000_1990 = calculateLEI(new0014, existing1990, dist=300)
                    #Attache LEI results to dataframe    
                    new0014['LEI'] = pd.DataFrame(lei_20142000_1990)
                    new0014.to_file(os.path.join(outFolder, "LEI_20142000_1990.shp"))
            except: 
                pass
            




