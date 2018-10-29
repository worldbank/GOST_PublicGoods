################################################################################
# Prepare and run Urban Analysis Data
# Benjamin P. Stewart
# Purpose: Prepare all the data necessary for urban analysis
#
# To Do: Ensure raster layers are extracted in map preparation process
################################################################################

import os, sys, inspect, argparse, logging, shutil

import geopandas as gpd
import pandas as pd

try:
    import rasterio
    from gbdxtools import Interface
    from shapely.geometry import Point, LineString
    from shapely.wkt import loads
    from shapely.geometry import box
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
    from GOSTRocks import rasterMisc
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

class urbanAnalysis_arcpy:
    '''certaing functions require older ArcPy analysis that needs to run separately'''
    def __init__(self, **kwargs):
        #set input variables
        self.inShapefile = kwargs['inShape']
        if 'urbanFolder' in kwargs.keys():
            self.urbanFolder = kwargs['urbanFolder']
        else:
            self.urbanFolder = os.path.dirname(self.inShapefile)
        self.gridFile = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_GRID.shp"))        
        self.gridcsv = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_GRID.csv"))
        self.lcFile = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_lcSummary.csv"))
        self.imageryFile = os.path.join(self.urbanFolder, os.path.basename(self.gridFile).replace(".shp", "_imagery_search.csv"))
        self.osmSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_OSM_Summary.csv"))
        self.ghslSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_ghsl_Summary.csv"))
        
        #Define files for market access analysis
        self.healthFacilities = os.path.join(self.urbanFolder, "OSM_Health.csv")
        self.healthMatrix = os.path.join(self.urbanFolder, "OSM_GRID_Health_Matrix.csv")
        self.healthMA = os.path.join(self.urbanFolder, "OSM_GRID_Health_MarketAccess.csv")
        self.eduFacilities = os.path.join(self.urbanFolder, "OSM_Education.csv")
        self.eduMatrix = os.path.join(self.urbanFolder, "OSM_GRID_Education_Matrix.csv")
        self.eduMA = os.path.join(self.urbanFolder, "OSM_GRID_Education_MarketAccess.csv")
        
        #Define files for mapping
        self.mappingFolder = os.path.join(self.urbanFolder, "MappingData")
        self.mapsFolder = os.path.join(self.urbanFolder, "Maps")
        self.healthMapping = os.path.join(self.mappingFolder, "HealthPoints.shp")
        self.eduMapping = os.path.join(self.mappingFolder, "EducationPoints.shp")
        self.gridMapping = os.path.join(self.mappingFolder, "grid_file.shp")
        self.fluvialMap = os.path.join(self.mappingFolder, "FU1_Combined.tif")
        self.pluvialMap = os.path.join(self.mappingFolder, "PU1_Combined.tif")
        self.combinedFlood = os.path.join(self.mappingFolder, "FU_PU_1_Combined.tif")
        self.mapFile = os.path.join(self.mappingFolder, "City_Scan.mxd")
        
    def prepareMaps(self):
        ''' Add input data to the map and symbolize
        '''
        studyArea = {'Study Area': self.inShapefile,
                    'Education Points': self.eduMapping,
                    'Education MA': self.gridMapping,
                    'Health Points': self.healthMapping,
                    'Health MA': self.gridMapping,
                    'OSMLR Level 1': self.gridMapping,
                    'OSMLR Level 2': self.gridMapping,
                    'OSMLR Level 3': self.gridMapping,
                    'OSMLR Level 4': self.gridMapping,
                    'GHSL Total Built': self.gridMapping,
                    'NDVI': self.gridMapping,
                    'BAI': self.gridMapping,
                    'Elevation': self.gridMapping,
                    'ESA AfriCover': os.path.join(self.mappingFolder, "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif.tif"),
                    'GHSL': os.path.join(self.mappingFolder, "GHSL.tif"),
                    'Fluvial Flood': self.fluvialMap,
                    'Pluvial Flood': self.pluvialMap,
                    'Combined Flood': self.combinedFlood
                    }
        #Open the mapFile
        #if not os.path.exists(self.mapFile):
        baseMap = os.path.join(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])), "City_Scan.mxd")
        shutil.copyfile(baseMap, self.mapFile)
            
        mxd = arcpy.mapping.MapDocument(self.mapFile)
        mxd.relativePaths = False
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, '', df)        
        for lyr in layers:
            try:
                newPath = studyArea[lyr.name]
                try:
                    lyr.replaceDataSource(os.path.dirname(newPath), "SHAPEFILE_WORKSPACE", os.path.basename(newPath).replace(".shp", ""))
                except:
                    lyr.replaceDataSource(os.path.dirname(newPath), "RASTER_WORKSPACE")
                
            except:
                print("Mapping %s: could not update data source" % lyr.name)
            if lyr.name == "Study Area":
                arcpy.SelectLayerByAttribute_management(lyr, "NEW_SELECTION", '"FID" = 0')
                df.extent = lyr.getSelectedExtent()
                arcpy.SelectLayerByAttribute_management(lyr, "CLEAR_SELECTION", '')                                       
        #Search for GHSL tiles and add to map        
        
        mxd.save()
    
    def generateMaps(self):
        ''' Create output mapsTODO: Add and symbolize layers
            2. If map doesn't exist, copy it from somewhere
        '''
        mxd = arcpy.mapping.MapDocument(self.mapFile)
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "Education.png"),  visibility = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "Hospital.png"),   visibility = [1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "OSM_Level1.png"), visibility = [1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "OSM_Level2.png"), visibility = [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "OSM_Level3.png"), visibility = [1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "OSM_Level4.png"), visibility = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "GHSL_Total.png"), visibility = [1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "NDVI.png"),       visibility = [1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "BAI.png"),        visibility = [1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1])
        createMapFromMXD(mxd, os.path.join(self.mapsFolder, "Elevation.png"),  visibility = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1])
    
    
class urbanAnalysis:
    def __init__(self, **kwargs): #*args
        self.wgs84 = {'init': u'epsg:4326'}
        self.webMercator = {'init': u'epsg:3857'}
        #set input variables
        self.inShapefile = kwargs['inShape']
        self.name = os.path.basename(self.inShapefile).replace(".shp", "")
        self.s3bucket = 'bps/cityAnalysis/%s' % self.name        
        self.inputAOI = gpd.read_file(self.inShapefile)
        self.inputAOI_WGS84 = self.inputAOI.to_crs(self.wgs84)        
        self.inputAOI_WM = self.inputAOI.to_crs(self.webMercator)        
        if 'urbanFolder' in kwargs.keys():
            self.urbanFolder = kwargs['urbanFolder']
        else:
            self.urbanFolder = os.path.dirname(self.inShapefile)            
        if 'osrmHeader' in kwargs.keys():
            self.osrmHeader = kwargs['osrmHeader']
        else:
            self.osrmHeader = 'http://router.project-osrm.org/table/v1/driving/'
        self.gridFile = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_GRID.shp"))        
        self.gridcsv = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_GRID.csv"))
        self.lcFile = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_lcSummary.csv"))
        self.imageryFile = os.path.join(self.urbanFolder, os.path.basename(self.gridFile).replace(".shp", "_imagery_search.csv"))
        self.osmSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_OSM_Summary.csv"))
        self.ghslSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_ghsl_Summary.csv"))
        self.ndviSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_NDVI_Summary.csv"))    
        self.baiSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_BAI_Summary.csv"))    
        self.srtmSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_SRTM_Summary.csv"))    
        self.SSBNSummary = os.path.join(self.urbanFolder, os.path.basename(self.inShapefile).replace(".shp", "_SSBN_Summary.csv"))    
        self.SSBNReturn = self.SSBNSummary.replace(".csv", "_RETURNDate.csv")
        
        #Define files for market access analysis
        self.healthFacilities = os.path.join(self.urbanFolder, "OSM_Health.csv")
        self.healthMatrix = os.path.join(self.urbanFolder, "OSM_GRID_Health_Matrix.csv")
        self.healthMA = os.path.join(self.urbanFolder, "OSM_GRID_Health_MarketAccess.csv")
        self.eduFacilities = os.path.join(self.urbanFolder, "OSM_Education.csv")
        self.eduMatrix = os.path.join(self.urbanFolder, "OSM_GRID_Education_Matrix.csv")
        self.eduMA = os.path.join(self.urbanFolder, "OSM_GRID_Education_MarketAccess.csv")
        
        #Define files for mapping
        self.mappingFolder = os.path.join(self.urbanFolder, "MappingData")
        self.mapsFolder = os.path.join(self.urbanFolder, "Maps")
        self.healthMapping = os.path.join(self.mappingFolder, "HealthPoints.shp")
        self.eduMapping = os.path.join(self.mappingFolder, "EducationPoints.shp")
        self.OSMRoads = os.path.join(self.mappingFolder, "OSM_Roads.shp")
        self.gridMapping = os.path.join(self.mappingFolder, "grid_file.shp")
        self.ghsl_raw = os.path.join(self.mappingFolder, "GHSL.tif")
        
        for cFolder in [self.mappingFolder, self.mapsFolder]:
            if not os.path.exists(cFolder):
                os.mkdir(cFolder)
                
        #Create grid file
        if not os.path.exists(self.gridFile):
            gridSize = 1000
            if 'gridSize' in kwargs.keys():
                gSize = kwargs['gridSize']
                #if len(gSize) > 1:
                #    gSize = gSize[0]
                gridSize = int(gSize)
            self.createGrid(self.gridFile, gridDist = gridSize)        

        #Get grid CRS
        self.gridGPD = gpd.read_file(self.gridFile)
        self.crs = self.gridGPD.crs              
    
    def checkFileStatus(self, inFile, name):
        if os.path.exists(inFile):
            print("%s has been PROCESSED" % name)
        else:
            print("***** %s has not been PROCESSED" % name)            
    
    def checkStatus(self):
        #Check for outputFolders
        self.checkFileStatus(self.lcFile, "LANDCOVER")
        self.checkFileStatus(self.imageryFile, "IMAGERY SEARCH")
        self.checkFileStatus(self.osmSummary, "OSM SUMMARY")
        #Check column names for grids to evaluate GHSL summary
        if "totalNB" in self.gridGPD.columns:
            print("%s has been PROCESSED" % "GHSL")
        else:
            print("***** %s has not been PROCESSED" % "GHSL")
    
    def summarizeGHSL(self):
        if not os.path.exists(self.ghslSummary):
            logging.info("Running GHSL Summary")            
            urbanParams = misc.getUrbanParams()
            ghslVRT = urbanParams["ghslVRT"]
            ghslRes = rasterMisc.zonalStats(self.gridFile, ghslVRT, reProj = True, 
                        verbose=True , rastType='C', unqVals=[0,1,2,3,4,5,6])
            ghslRes = pd.DataFrame(ghslRes, columns=['NoData','Water','NotBuilt','b2014','b2000','b1990','b1975'])
            #Create total built metrics
            ghslRes['TotalBuilt'] = ghslRes['b2014'] + ghslRes['b2000'] + ghslRes['b1990'] + ghslRes['b1975']
            ghslRes.to_csv(self.ghslSummary)  
        else:
            logging.info("GHSL Summary already exists")            
    
    def clipGHSL(self):
        urbanParams = misc.getUrbanParams()
        ghslVRT = urbanParams["ghslVRT"]
        rasterMisc.clipRaster(rasterio.open(ghslVRT), self.gridGPD, self.ghsl_raw)
    
    def summarizeLULC(self, lulcFile, outFile):
        ''' Run zonal stats against input lulcFile; such as the GBDx 
        '''
        lulcRes = rasterMisc.zonalStats(self.gridFile, lulcFile, reProj=True, verbose=True, rastType='C', 
                                        unqVals = [0,1,2,3,4,5,6])
        lulcRes = pd.DataFrame(lulcRes, columns=["NoData","Veg","H2O","Bare","Cloud","Shdw","Built"])
        lulcRes.to_csv(lulcFile)
    
    def summarizeSSBN(self, ssbnFolder, imageType="tif"):
        ''' Run zonal statistics on all files in the input SSBN folder
        '''
        #Get reference to the appropriate flood files
        if ssbnFolder == '':
            raise ValueError("No SSBN Folder defined")
        if not os.path.exists(self.SSBNSummary):
            #Get a reference to SSBN flood tif files
            allTifs = []
            for root, dirs, files in os.walk(ssbnFolder):
                for f in files:
                    if f[-3:] == imageType:
                        allTifs.append(os.path.join(root, f))
            #Run zonal stats on all SSBN tif files
            for cTif in allTifs:
                name = os.path.basename(cTif).replace(imageType, "")
                logging.info("Processing SSBN file %s" % name)
                ##Get the flood type from the raster name
                columnNames = ["%s-%s" % (name, x) for x in ["SUM","MIN","MAX","MEAN"]]
                curRes = rasterMisc.zonalStats(self.gridFile, cTif, rastType='N', reProj=True, verbose=False)
                curPD = pd.DataFrame(curRes, columns=columnNames)
                #Drop the Min and Sum columns
                curPD = curPD.drop(curPD.columns[[0,1]], axis=1)
                try:
                    final = final.join(curPD)
                except:
                    final = curPD
            final.to_csv(self.SSBNSummary)                   
        else:
            final = pd.read_csv(self.SSBNSummary)        
        #Create Flood derived statistics
        #Which cells have both fluvuial and pluvial risk
        returnRates = [1000,500,250,200,100,75,50,20,10,5]
        curStat = "MAX"                

        def calculateMinReturn(x, minVal=0, maxVal=100):
            ''' For each row, identify the column where the first non-0 value appears and return that year '''
            fVals = x[(x > minVal) & (x < maxVal)]
            try:
                return(min([int(x.split("-")[2]) for x in fVals.index]))
            except:
                return 0

        #For each cell, determine the return rate at which it floods
        #Select appropriate columns for each of U and D
        curStat = "MAX"  
        curD = final[[col for col in final.columns if curStat in col]]
        cFlood = "FD"
        allRes = []
        for cFlood in ["FD", "FU", "PD", "PU", "UD", "UU"]:
            curD = final[[col for col in final.columns if cFlood in col]]
            allRes.append(curD.apply(calculateMinReturn, axis=1))

        xx = pd.DataFrame(allRes).transpose()
        xx.columns = ["FD", "FU", "PD", "PU", "UD", "UU"]
        xx.to_csv(self.SSBNReturn)               
    
    def createGrid(self, outGrid, gridDist=1000):
        ''' The urban analysis is run on a grid within the define shapefile AOI
            [Variables]
            outGrid - output grid shapefile
            gridDist - size (in m) of output grid
        '''
        logging.info("Creating input grid")            
        fGrid = self.inputAOI_WM
        if fGrid.crs == {'init': u'epsg:4326'}:
            logging.warn("Input Shapefile is in WGS84, which could mess up the grid if not careful")
        b = fGrid.unary_union.bounds
        curCrs = int(self.inputAOI_WM.crs['init'].split(":")[1])
        misc.createFishnet(self.gridFile,b[0], b[2], b[1], b[3], gridDist, gridDist, crsNum=curCrs)
        #Eliminate grids that do not intersect study area
        curGPD = gpd.read_file(self.gridFile)
        curGPD = curGPD[curGPD.intersects(self.inputAOI_WM.unary_union)]
        curGPD['FID'] = range(0, curGPD.shape[0])
        curGPD.to_file(self.gridFile)
    
    def summarizeOSM(self, verbose=True):
        ''' There are a number of metrics to be extracted from OSM within the urban grid
            1. Density/Sum of length of roads
            2. Density/Sum of railroads
            3. Market access (to airports?)            
        if not os.path.exists(self.osmSummary):
            xx = osmMisc.summarizeOSM(self.gridGPD, verbose=verbose)
            xx.to_csv(self.osmSummary)
        else:
            logging.info("OSM Summary already exists")            
        '''
        if not os.path.exists(self.OSMRoads):
            logging.info("Creating OSM Road shapefile")            
            xx = osmMisc.summarizeOSM(self.gridGPD, verbose=verbose)
            xx.to_file(self.OSMRoads, roadsOnly=False)
            
    
    def summarizeElevation(self):
        ''' Summarize SRTM elevation within grid using Truele Earth Engine        
        '''
        #Summarize BAI
        if not os.path.exists(self.srtmSummary):  
            logging.info("Running SRTM Elevation statistics")
            #Project the data to Web Mercator, if it isn't
            baiGrid = self.gridFile
            if self.crs != self.wgs84:
                curD = self.gridGPD.to_crs(self.wgs84)
                baiGrid = self.gridFile.replace(".shp", "_WGS84.shp")
                curD.to_file(baiGrid)            
            #Run google earth engine with the current shapefile   
            baiImages = gee.getElevimages()
            gee.zonalStats(baiGrid, baiImages, self.srtmSummary)
        else:
            logging.info("Elevation statistics already created")
            
    def summarizeBAI(self):
        ''' Summarize Built Area Index within grid using Google Earth Engine        
        '''
        #Summarize BAI
        if not os.path.exists(self.baiSummary):  
            logging.info("Running BAI statistics")
            #Project the data to Web Mercator, if it isn't
            baiGrid = self.gridFile
            if self.crs != self.wgs84:
                curD = self.gridGPD.to_crs(self.wgs84)
                baiGrid = self.gridFile.replace(".shp", "_WGS84.shp")
                curD.to_file(baiGrid)            
            #Run google earth engine with the current shapefile   
            baiImages = gee.getBAIimages(largeTimeSteps=True)
            gee.zonalStats(baiGrid, baiImages, self.baiSummary)
        else:
            logging.info("BAI statistics already created")
    
    def summarizeNDVI(self):
        ''' Summarize NDVI within grid using Google Earth Engine        
        '''
        #Summarize NDVI
        if not os.path.exists(self.ndviSummary):  
            logging.info("Running NDVI statistics")
            #Project the data to Web Mercator, if it isn't
            ndvi_grid = self.gridFile
            if self.crs != self.wgs84:
                curD = self.gridGPD.to_crs(self.wgs84)
                ndvi_grid = self.gridFile.replace(".shp", "_WGS84.shp")
                curD.to_file(ndvi_grid)            
            #Run google earth engine with the current shapefile   
            ndviImages = gee.getNDVIimages(largeTimeSteps=True)
            gee.zonalStats(ndvi_grid, ndviImages, self.ndviSummary)
        else:
            logging.info("NDVI statistics already created")
            
    def searchForImagery(self):
        '''search for imagery within the defined AOI
        '''
        if not os.path.exists(self.imageryFile):
            logging.info("Saerching for imagery within AOI")            
            gbdx = Interface()
            curWKT = self.inputAOI_WGS84['geometry'][0]
            curRes = imagery_search.searchForImages(gbdx, curWKT, "C:/Temp", "Balikpanana", 
                                cutoff_date='1-Jan-15', optimal_date='17-July-18',
                                cutoff_cloud_cover = 100, cutoff_nadir = 90)
            curRes.to_csv(self.imageryFile)
        else:
            logging.info("Imagery file already exists")            

    def summarizeLandcover(self, lcFile, unqVals = [1,2,3,4,5,6,7,8,10,200]):
        '''run categorical zonal stats on lcFile
        ###TODO: calculate zonal values from the input LCFile
        '''
        if not os.path.exists(self.lcFile):
            logging.info("Running categorical summary of landcover file")            
            xx = rasterMisc.zonalStats(self.gridFile, lcFile, rastType='C', reProj=True, verbose=True,unqVals=unqVals)
            outX = pd.DataFrame(xx, columns=unqVals)
            outX.to_csv(self.lcFile)
        else:
            logging.info("Landcover summary already exists") 
            
    def MarketAccess(self):
        ''' Run market analysis on OSM POIs
        '''
        #Check if POIs already exist - create them if they don't
        if not os.path.exists(self.healthFacilities):            
            current = OSMNX_POIs.AmenityObject('Health', self.inputAOI_WGS84.unary_union, ['clinic','pharmacy','hospital','health'], self.urbanFolder)
            temp = current.GenerateOSMPOIs()
            temp = current.RemoveDupes(0.005, self.inputAOI_WGS84.crs)
            temp = current.prepForMA()
            temp['TotalBuilt'] = 1
            temp.to_csv(self.healthFacilities, encoding='utf-8')
        if not os.path.exists(self.eduFacilities):
            current = OSMNX_POIs.AmenityObject('Education', self.inputAOI_WGS84.unary_union, ['school','university','secondary school', 'kindergarten', 'college'], self.urbanFolder)
            temp = current.GenerateOSMPOIs()
            temp = current.RemoveDupes(0.005, self.inputAOI_WGS84.crs)
            temp = current.prepForMA()
            temp['TotalBuilt'] = 1
            temp.to_csv(self.eduFacilities, encoding='utf-8')            
        if not os.path.exists(self.gridcsv):
            #Convert the grid data into points
            ptGPD = self.gridGPD
            ptGPD = ptGPD.to_crs(self.wgs84)
            ptGPD.geometry = ptGPD.centroid
            ptGPD['Lon'] = [x.x for x in ptGPD.geometry]
            ptGPD['Lat'] = [x.y for x in ptGPD.geometry]
            #Attach GHSL information to be used as weight
            if not os.path.exists(self.ghslSummary):
                self.summarizeGHSL()                
            ghslPD = pd.read_csv(self.ghslSummary)
            ptGPD = pd.concat([ptGPD, ghslPD], axis=1)
            ptGPD['mID'] = range(0, ptGPD.shape[0])
            ptGPD.to_csv(self.gridcsv, encoding='utf-8')
        
        #Run OD Matrix from grid to each of the facilities
        if not os.path.exists(self.healthMatrix):
            healthMatrix = OD.CreateODMatrix(self.gridcsv,self.healthFacilities, Pop='TotalBuilt', UID='mID', 
            sleepTime=0, osrmHeader=self.osrmHeader)
            healthMatrix.to_csv(self.healthMatrix)
        else:
            healthMatrix = pd.read_csv(self.healthMatrix)
        if not os.path.exists(self.healthMA):
            healthMA = OD.MarketAccess(healthMatrix)
            healthMA.to_csv(self.healthMA)
            
        if not os.path.exists(self.eduMatrix):
            eduMatrix = OD.CreateODMatrix(self.gridcsv, self.eduFacilities, Pop='TotalBuilt', UID='mID', 
                sleepTime=0, osrmHeader=self.osrmHeader)
            eduMatrix.to_csv(self.eduMatrix)
        else:
            eduMatrix = pd.read_csv(self.eduMatrix)
        if not os.path.exists(self.eduMA):
            eduMA = OD.MarketAccess(eduMatrix)
            eduMA.to_csv(self.eduMA)
            
    
    def summarizeNumberOfCars(self, catalogID, outFolder):
        ''' Run car counting from GBDx against imageyr for the defined area       
        '''
        
    def runSPFEAS(self, catIDs, sensorList):
        ''' Summarize spfeas through zonal statistics
        '''
        gbdx = Interface()
        curTasks = gbdxTasks.GOSTTasks(gbdx)
        gbdxUrl = gbdxURL_misc.gbdxURL(gbdx)
        allWorkflows = []
        curGeom = str(self.gridGPD.to_crs(self.wgs84).unary_union)
        for catID_idx in range(0, len(catIDs)):
            catID = catIDs[catID_idx]
            sensor = sensorList[catID_idx]
            data = gbdx.catalog.get_data_location(catID)
            if data == 'not_delivered':
                #order the image if it is ont ordered yet
                #orderID = gbdx.ordering.order(catID)
                #curStatus = gbdx.ordering.status(orderID)
                #print curStatus
                print ("Waiting for image ordering - %s" % catID)
            else:
                #Run spfeas on analysis
                print ("Executing spfeas and car counting on %s" % catID)
                allWorkflows.append(curTasks.createWorkflow(catID, curGeom, 
                    sensor, "%s/%s" % (self.s3bucket, catID),
                    spfeasParams={"triggers":'dmp fourier gabor grad hog lac lsr mean pantex orb saliency sfs', 
                                  "scales":'8 16 32', "block":'4', "gdal_cache":64, "section_size":2000}, 
                    runCarFinder = 0, runSpfeas = 1, downloadImages = 1,
                    aopPan=False, aopDra=False, aopAcomp=True, aopBands='MS'))
        for a in allWorkflows:
            a.execute()
        xx = gbdxUrl.monitorWorkflows(sleepTime=300)
    
    def prepMappingData(self, floodFolder=''):
        '''convert all the data to shapefiles for use in mapping
           1. Health / Education points
           2. Health / Education access / routes
           3. OSM Summary
           4. LC Summary / Landcover
           5. GHSL Summary / GHSL
           
           Clip input raster datasets           
        '''
        
        def createPoints(inFile, outFile):
            #Map OSM point files
            if not os.path.exists(outFile):
                inDF = pd.read_csv(inFile)
                geoms = [Point(xy) for xy in zip(inDF.Lon, inDF.Lat)]
                inDF = inDF.drop(['Lat', 'Lon'], axis=1)
                inGDF = gpd.GeoDataFrame(inDF, geometry=geoms, crs=self.wgs84)
                inGDF.to_file(outFile)
        
        def mergeFile(inD, inFile, prefix):
            try:
                mData = pd.read_csv(inFile)
                mData.columns = [prefix + s.replace(".", "") for s in mData.columns]
                inD = inD.merge(mData, how='left', left_on="FID", right_on=mData.columns[0])
            except:
                logging.info("could not process %s" % inFile)
            return(inD)
                    
        createPoints(self.healthFacilities, self.healthMapping)
        createPoints(self.eduFacilities, self.eduMapping)
        if not os.path.exists(self.gridMapping):
            #A number of attributes need to be attached to the grid file
            inGrid = self.gridGPD
            #Merge the market access files, the osm summary, the ghsl summary
            inGrid = mergeFile(inGrid, self.healthMA, "H_")
            inGrid = mergeFile(inGrid, self.eduMA, "E_")
            inGrid = mergeFile(inGrid, self.osmSummary, "OSM_")
            inGrid = mergeFile(inGrid, self.ghslSummary, "GHSL_")
            inGrid = mergeFile(inGrid, self.lcFile, "LC_")
            inGrid = mergeFile(inGrid, self.ndviSummary, "NDVI_")
            inGrid = mergeFile(inGrid, self.baiSummary, "BAI_")
            inGrid = mergeFile(inGrid, self.srtmSummary, "ELEV_")
            inGrid = mergeFile(inGrid, self.SSBNReturn, "SSBN_")
            inGrid.to_file(self.gridMapping)
        
        #Clip input raster datasets
        urbanParams = misc.getUrbanParams()
        #Get reference to combined flood data
        curRasters = [urbanParams['afriCover20'],urbanParams['ghspop15'],urbanParams['ghslVRT'],
                        urbanParams['worldPopAfrica']]
        if floodFolder != '':
            cFloodFolder = os.path.join(floodFolder, "CombinedFloodData")
            floodRasters = os.listdir(cFloodFolder)
            for r in floodRasters:
                try:
                    curR = rasterio.open(os.path.join(cFloodFolder, r), 'r')
                    floodAOI = self.inputAOI.to_crs(curR.crs)
                    intersects = box(curR.bounds.left,curR.bounds.bottom,curR.bounds.right,curR.bounds.top).intersects(floodAOI.unary_union)
                    if intersects:                    
                        curRasters.append(os.path.join(cFloodFolder, r))                                
                except:
                    pass
        for inRaster in curRasters:
            outRaster = os.path.join(self.mappingFolder, os.path.basename(inRaster).replace(".vrt", ".tif"))  
            if not os.path.exists(outRaster):
                rasterMisc.clipRaster(rasterio.open(inRaster), self.gridGPD, outRaster) 

if __name__ == "__main__":
    exampleText =  '''
    #Run OSRM with custom OSRM link 
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\Indonesia_GBDx\BalikPapan_AOI.shp -osm -ma -osrm http://router.project-osrm.org/table/v1/driving/
    
    #Run SPFEAS analysis
    python UrbanAnalysis.py -i Q:\AFRICA\COD\Projects\DRC_Census_Franck\Goma_WM\Goma_WM.shp -spfeas 105001000AC23A00 10400100091AC900 -sensor GEOEYE01 WORLDVIEW03_VNIR
    
    #Run Landcover Analysis
    python UrbanAnalysis.py -i Q:\AFRICA\COD\Projects\DRC_Census_Franck\Goma_WM\Goma_WM.shp -lcFile Q:\AFRICA\LandCover\Africover_2016_20m\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif    

    #Run Google Earth Engine analysis
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\ImageryDownload\Mali_Keith\Ansongo_AOI.shp -gdist 250 -ndvi -bai   
    
    #Run all basic analysis
    python UrbanAnalysis.py -i Q:\AFRICA\STP\Projects\CEM\SaoTome.shp -bai -ndvi -elev -ghsl -lcFile Q:\AFRICA\LandCover\Africover_2016_20m\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif -osm -ma

    #run all analyses with custom grid size
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\CityScan\BalikPapan_AOI.shp -ghsl -lcFile Q:\GLOBAL\LCVR\Globcover\2015\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif -osm -ma -gdist 250

    #Run SSBN analysis
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\CityScan\Data\CityExtents\Conotou_AOI.shp -ssbn Q:\GLOBAL\HYDRO\SSBN_Flooding\benin
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\Indonesia_GBDx\BalikPapan_AOI.shp -ssbn Q:\GLOBAL\HYDRO\SSBN_Flooding\indonesia
    python UrbanAnalysis.py -i Q:\WORKINGPROJECTS\CityScan\Data\CityExtents\Douala_AOI.shp -ssbn Q:\GLOBAL\HYDRO\SSBN_Flooding\cameroon
    
    ###RECENT Testing
    python UrbanAnalysis.py -i Q:\AFRICA\COD\Projects\DRC_Census_Franck\Kinshasai_WM\Kinshasai.shp -gdist 100 -ghsl 
    python UrbanAnalysis.py -i Q:\AFRICA\COD\Projects\DRC_Census_Franck\Kinshasai_WM\Kinshasai.shp -gdist 100 -lcFile Q:\AFRICA\LandCover\Africover_2016_20m\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif 
    python UrbanAnalysis.py -i Q:\AFRICA\COD\Projects\DRC_Census_Franck\Kinshasai_WM\Kinshasai.shp -ma -osrm http://w1es235:5000/table/v1/driving/

    '''
    parser = argparse.ArgumentParser(description="Generate urban metrics based on a defined grid")
    
    parser.add_argument('-i',       dest="INPUT_AOI", action='store', help="Shapefile AOI")
    parser.add_argument('-o',       dest="OUT_FOLDER", action='store', help="outputFolder")
    parser.add_argument('-lcFile',  dest='LCFILE', action='store', help="Landcover file to process - defaults to Globcover")
    parser.add_argument('-spfeas',  dest='SPFEAS', action='store', help="Catalogue IDs for running spfeas", nargs='+')
    parser.add_argument('-sensor',  dest='SENSOR', action='store', help="Sensor for processing spfeas, default is WORLDVIEW03_VNIR", default="WORLDVIEW03_VNIR", nargs='+')
    
    parser.add_argument('-gdist',   dest='GDIST', action='store', help="spacing for grid overlaid over dataset", default=1000)
    parser.add_argument('-arcpy',   dest='ARCPY', action='store_true', help="Create output maps **Note** this needs to be run from a python instance with access to ArcPy")
    parser.add_argument('-s',       dest='STATUS', action='store_true', help="Check on status of processing")
    parser.add_argument('-map',     dest='MAP', action='store_true', help="Create data for mapping")
    parser.add_argument('-ghsl',    dest='GHSL', action='store_true', help="Summarize GHSL within grid")
    parser.add_argument('-ndvi',    dest='NDVI', action='store_true', help="Summarize NDVI within grid")
    parser.add_argument('-bai',     dest='BAI', action='store_true', help="Summarize BAI within grid")
    parser.add_argument('-elev',    dest='ELEV', action='store_true', help="Summarize SRTM elevation within grid")
    parser.add_argument('-osm',     dest='OSM', action='store_true', help="Summarize OSM Road Lengths")
    parser.add_argument('-osrm',    dest='OSRMHEADER', action='store', help="OSRM base path for running market access")
    parser.add_argument('-ma',      dest='MA', action='store_true', help="Run Market Analysis")
    parser.add_argument('-ssbn',    dest='SSBNFolder', action='store', help="Folder containing the SSBN results")
    
    args = parser.parse_args()
    #Set logging information
    logging.basicConfig(format='%(asctime)s:%(levelname)s:  %(message)s', level=logging.INFO)
    
    inFile = args.INPUT_AOI
    try:
        outFolder = inFile.replace(".shp", "")
        if args.OUT_FOLDER:
            outFolder = args.OUT_FOLDER
        if not os.path.exists(outFolder):
            os.mkdir(outFolder)
    except:
        pass
    
    
    if args.ARCPY:
        urbanAnalysis_arcpy(inShape=inFile, urbanFolder=outFolder)
    else:        
        ua = urbanAnalysis(inShape=inFile, urbanFolder=outFolder, gridSize=int(args.GDIST))               
        
        if args.STATUS:
            ua.checkStatus()
        if args.MAP:
            ua.prepMappingData()
        if not os.path.exists(ua.imageryFile):
            ua.searchForImagery()             
        if args.SSBNFolder:
            ua.summarizeSSBN(args.SSBNFolder) 
        if args.NDVI:
            ua.summarizeNDVI()        
        if args.BAI:
            ua.summarizeBAI()
        if args.ELEV:
            ua.summarizeElevation()
        if args.LCFILE:
            ua.summarizeLandcover(args.LCFILE)
        if args.OSM:
            ua.summarizeOSM()
        if args.GHSL:
            ua.summarizeGHSL()
        if args.MA:
            ua.MarketAccess()
        if args.SPFEAS:
            sensor = args.SENSOR
            if len(sensor) != len(args.SPFEAS):
                sensor = [sensor] * len(args.SPFEAS)
            ua.runSPFEAS(args.SPFEAS, sensor)

