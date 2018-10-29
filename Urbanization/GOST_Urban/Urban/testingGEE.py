###TESTING GEE
import os, sys, inspect, argparse, logging, shutil, re, geojson
import ee
ee.Initialize() #Initialize the ee object, but this requires work the first time through 
                # https://developers.google.com/earth-engine/python_install

import geopandas as gpd
import pandas as pd

try:
    from gbdxtools import Interface
    from shapely.geometry import Point, LineString
    from shapely.wkt import loads
except:
    pass
#Import a number of existing GOST functions
GOSTRocks_folder = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))))
GOST_GBDx_folder = os.path.join(os.path.dirname(GOSTRocks_folder), "GOST_GBDx")
GOST_Public_Goods = os.path.join(os.path.dirname(GOSTRocks_folder), "GOST_PublicGoods")

if GOSTRocks_folder not in sys.path:
    sys.path.insert(0, GOSTRocks_folder)

if GOST_GBDx_folder not in sys.path:
    sys.path.insert(0, GOST_GBDx_folder)

if GOST_Public_Goods not in sys.path:
    sys.path.insert(0, GOST_Public_Goods)

from GOSTRocks import misc
try:
    from Market_Access import OD
    from Market_Access import OSMNX_POIs
    from GOSTRocks import osmMisc
    from GOSTRocks import GEE_zonalstats as gee
    from GOST_GBDx_Tools import gbdxTasks
    from GOST_GBDx_Tools import gbdxURL_misc
    from GOST_GBDx_Tools import imagery_search
except:
    pass
try:
    import arcpy
    from GOSTRocks.Urban.SummarizeGHSL import *
    from GOSTRocks.arcpyMisc import createMapFromMXD
    ARCPY_LOADED = True
except:
    pass

logging.basicConfig(level=logging.INFO)    
inputShapefile = r"Q:\WORKINGPROJECTS\Indonesia_GBDx\BalikPapan_AOI\DEBUG\BalikPapan_AOI_GRID_WGS84.shp"
baiImages = []
baiImages.append(ee.Image("LANDSAT/LC8_L1T_ANNUAL_BAI/%s" % 2017))
baiImages.append(ee.Image("LANDSAT/LE7_L1T_ANNUAL_BAI/%s" % 2010))
baiImages.append(ee.Image("LANDSAT/LE7_L1T_ANNUAL_BAI/%s" % 2000))
        
xx = gee.zonalStats(inputShapefile, baiImages, inputShapefile.replace(".shp", ".csv"))
print xx


'''
inputImages = baiImages
inD = gpd.read_file(inputShapefile)
allFeatures = []
idx = 0
subIdx = 0
vCount = 0
scale=50
#Loop through the features in the shapefile
inShapes = inD['geometry']
shp = inShapes[1]
for shp in inShapes[0:3]:
    idx = idx + 1  
    subIdx = subIdx + 1
    shpJSON = geojson.Feature(geometry=shp, properties={})
    outputAttributes = {'index': subIdx}       
    for ndviImage in inputImages:
        bandName = str(ndviImage.bandNames())
        bName = re.search('id(.*)', bandName)
        bName = bName.group(1).replace('"', '').replace(':', '').replace(' ', '')
        outputAttributes[bName] = ndviImage.reduceRegion(
            reducer = ee.Reducer.percentile([0,25,50,75,100], None, None, None, None),
            #reducer=ee.Reducer.sum(), 
            geometry=ee.Geometry.Polygon(shpJSON['geometry']['coordinates']), 
            scale=scale, 
            maxPixels=10e15, 
            bestEffort=True)   
    allFeatures.append(ee.Feature(None, outputAttributes))  


def fetchResults(allFeatures):
    allFeaturesCollection = ee.FeatureCollection(allFeatures)
    ndviDict = allFeaturesCollection.getInfo()
    ###YOU FUCKING IDIOT BEN!@&$( FUCTIONS HAVE TO FUCKING RETURN THINGS
    return ndviDict

def extractResults(xx):
    allRes = []
    valNames = []
    logVals = True
    for x in xx['features']:
        curRes = []
        if x == xx['features'][0]:
            #Write the header if this is the first feature
            allHead = []
            for xHead in x['properties'].iterkeys():
                allHead.append(xHead)
        for val in x['properties'].itervalues():
            if type(val) is dict:
                for y in val.values():
                    curRes.append(y)
                if logVals:
                    valNames = val.keys()
                    logVals=False
            else:
                curRes.append(str(val))
        allRes.append(curRes)
    allNames = ['IDX']
    for hIdx in allHead[1:]:
        for vName in valNames:
            allNames.append(hIdx + "_" + vName)
    return {"Header":allNames, "Results": allRes}

ndviDict = fetchResults(allFeatures)
del final
allRes = []
for curFeat in ndviDict['features']:
    curRes = pd.DataFrame(curFeat['properties'])
    allRes.append(curRes.unstack())

x = pd.DataFrame(allRes)
x.columns = ['__'.join(col) for col in x.columns.values]
    
    try:
        final.append(curRes, axis=0)
    except:
        final = curRes


pd.DataFrame(ndviDict['features'][0]['properties'])

n1 = extractResults(ndviDict)

n0 = ndviDict
      
    if idx > 100 or idx == inD.shape[0]:
        logging.info("Processing %s of %s" % (subIdx, len(inShapes)))
        vCount = vCount + 1
        curOut = outputFile.replace(".csv", "_%s.csv" % vCount)
        if not os.path.exists(curOut):
            try:
                ndviDict = fetchResults(allFeatures)
            except:
                logging.warning("***** Encountered GEE error, taking a break")
                time.sleep(60)
                ndviDict = fetchResults(allFeatures)
            if ndviDict and writeTemp:
                writeResults(ndviDict, curOut)
            else:
                pass
            pRes = extractResults(ndviDict)                
            if idx == subIdx:
                final = pRes['Results']
            else:
                for x in pRes['Results']:
                    final.append(x)                     
        #Reset for next iteration
        idx = 0
        allFeatures = []    
try:
    finalPD = pd.DataFrame(final, columns=pRes['Header'])
    #finalPD = finalPD.append(inD, axis=1)
    finalPD.to_csv(outputFile)
    #writeResults(final, outputFile)
    return finalPD
except:
    print final
'''