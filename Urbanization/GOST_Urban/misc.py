import os, sys, time, math, subprocess, inspect, json, glob

import geopandas as gpd
import shapely
import fiona 
import pyproj
import rasterio
import numpy
import ogr

from osgeo import osr
from collections import Counter
from math import ceil
from rtree import index
from shapely.geometry import shape
from affine import Affine
from rasterio.features import rasterize

def renameDir(dir):
    '''
    Renames all files in a directory
    Used to add .dat to all ENVI files, but not their headers
    '''
    for f in os.listdir(dir):
        title = os.path.basename(pathAndFilename)
        os.rename(pathAndFilename, os.path.join(dir, "%s.dat" % title))

def listFiles(inFolder, pattern):
    outFiles = []
    for f in os.listdir(inFolder):
        if f.endswith(pattern):
            outFiles.append(os.path.join(inFolder, f))
        elif os.path.isdir(os.path.join(inFolder, f)):
            outFiles = outFiles + listFiles(os.path.join(inFolder, f), pattern)
    return(outFiles)

'''prints the time along with the message'''
def tPrint(s):
    print("%s\t%s" % (time.strftime("%H:%M:%S"), s))

def round_to_1(x):
    return(round(x, -int(math.floor(math.log10(x)))))

# Create an interable range made with decimal point steps
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step    
'''
Get the index of a specific [val] within a list of histogram values
'''
def getHistIndex(hIdx, val, maxVal=2000):
    for h in range(0, len(hIdx)):
        curH = hIdx[h]
        if curH > val:
            return(lastH)        
        lastH = h
    return(len(hIdx) -1)

'''
Convert a list of values into a percent of total 
'''    
def getHistPer(inD):
    tSum = listSum(inD)
    for hIdx in range(0,len(inD)):
        inD[hIdx] = inD[hIdx] / tSum
    return(inD)
    

def generateVRT(inFiles, outVRT, 
        gdalVRTfunction=r"C:\Python27\ArcGIS10.3\Lib\site-packages\osgeo\gdalbuildvrt.exe", fileList="C:/Temp/VRTList.txt"):
    vrtFile = open(fileList, 'w')
    for f in inFiles:
        vrtFile.write("%s\n" % f)
    vrtFile.close()
    subprocess.call("%s -input_file_list %s %s" % (gdalVRTfunction, fileList, outVRT))
    os.remove(fileList)
    
def getUrbanParams():
    thisFolder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    inputParameters = "%s/Urban/urbanParameters.json" % thisFolder
    with open(inputParameters, 'r') as data_file:                 
        jsonParams = json.load(data_file)
    return jsonParams
    
def getParams():
    '''
    '''
    thisFolder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    inputParameters = "%s/importantParameters.json" % thisFolder    
    with open(inputParameters, 'r') as data_file:                 
        jsonParams = json.load(data_file)
    return jsonParams
    
def selectByIntersection(allAdmin, inD, exact=False):
    ''' Select features from allAdmin that intersect inD
    allAdmin [shapely geometry] - describes the feature of interest
    inD [geopandas dataframe] - features to be intersected
    exact (optional) [boolean] - perform exact geometry match, default is False, which
                                 will only calculate bounding box intersection
    RETURNS [geopandas dataframe] - subset of inD that intersects allAdmin
    EXAMPLE:
        inS = gpd.read_file(sourceFile)   
        inD = gpd.read_file(intersectingFile)        
        if inS.crs != inD.crs:
            inS = inS.to_crs(inD.crs)
        allAdmin = shapely.ops.cascaded_union(MultiPolygon(inS['geometry'].tolist()))
        selectByIntersection(allAdmin, inD, exact=False)
    TODO: Projection checks
    '''
    #Read the source boundaries into an rtree opbject
    idxAdmin = index.Index()    
    try:
        for i, shape in enumerate(allAdmin):
            idxAdmin.insert(i, shape.bounds)
    except:
        idxAdmin.insert(0, allAdmin.bounds)
    #Enumerate through the intersecting file
    intersectingList = []
    for i, row in inD.iterrows():
        interSection = list(idxAdmin.intersection(row['geometry'].bounds))
        if len(interSection):
            intersectingList.append(row)
            
    returnDF = gpd.GeoDataFrame(intersectingList) 

    if exact:
        selData = []
        for j, origRow in returnDF.iterrows():    
            xx = allAdmin.intersects(origRow['geometry'])
            print ("%s - %s" % (j, xx))
            if xx:
                selData.append(origRow)
            returnDF = gpd.GeoDataFrame(selData) 
    return returnDF

def tabulateUnq(unqResults, verbose=False, columnPrefix="c"):
    cnt = 0
    tResults = len(unqResults)
    allVals = []
    for r in unqResults:
        allVals.append(r[0].tolist())
    flattened = [val for sublist in allVals for val in sublist] 
    unq = np.unique(np.array(flattened)).tolist()
    #unqCols = ["c_%s" % xxx for xxx in unq]
    allRes = []
    for r in unqResults:               
        try:
            curRes = [0] * len(unq)
            for idx in range(0, len(r[0].tolist())):
                curRes[unq.index(r[0].tolist()[idx])] = r[1].tolist()[idx]
        except:
            print (r)
        allRes.append(curRes)
    return pd.DataFrame(allRes, columns=["%s_%s" % (columnPrefix, xxx) for xxx in unq])    
    
def createFishnet(outputGridfn,xmin,xmax,ymin,ymax,gridHeight,gridWidth,crsNum=4326):
    ''' Create a fishnet shapefile inside the defined coordinates 
    outputGridfn: output shapefile to hold the grid
    xmin,xmax,ymin,ymax: coordinates for the extent of the fishnet
    gridHeight,gridWidth: dimensions of the grid cells in the outpout fishnet
    '''    
    #Define projection
    testSR = osr.SpatialReference()
    res = testSR.ImportFromEPSG(crsNum)
    if res != 0:
        raise RuntimeError(repr(res) + ': could not import from EPSG')
    
    # convert sys.argv to float
    xmin = float(xmin)
    xmax = float(xmax)

    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        os.remove(outputGridfn)
    outDataSource = outDriver.CreateDataSource(outputGridfn)
    outLayer = outDataSource.CreateLayer(outputGridfn, testSR, geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Close DataSources
    outDataSource.Destroy()

'''
Get list of VIIRS files based on input filtering
'''  
def getVIIRSFiles(baseFolder=r"Q:\GLOBAL\NTL\VIIRS", years="all", months="all", tile=-1, retType="images"):
    if years == "all":
        years = [2012, 2013, 2014, 2015, 2016]
    if months == "all":
        months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    returnList = {}
    for yr in years:
        for mn in months:
            if tile == -1:
                curFolder = os.path.join(baseFolder, "*" + str(yr) + str(mn))
            else:
                curFolder = os.path.join(baseFolder, str(yr) + str(mn))
                curFolder = os.path.join(curFolder, "TILE%s" % tile)
                
            curFiles = glob.glob("%s\\*.tif" % curFolder)
            returnList[str(yr) + str(mn).zfill(2)] = curFiles
    
    if retType == "images":
        finalList = []
        for key, value in returnList.iteritems():
            if len(value) > 0:
                finalList.append(value[0])
        return(sorted(finalList))
    
    if retType == "clouds":
        finalList = []
        for key, value in returnList.iteritems():
            if len(value) > 0:
                finalList.append(value[1])
        return(sorted(finalList))
    return(returnList)
