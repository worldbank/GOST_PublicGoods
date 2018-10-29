#######################################################################
# Summarize Urban Clusters in Admin Districts
# Benjamin Stewart, October, 2017
# Purpose: compare admin boundaries and urban clusters
#   1. Input Layers
#       a. National admin boundaries
#       b. National population layer
#       c. National urban layer
#   2. Convert urban layers to vector
#       b. attach names from GRUMP
#   3. For each HDC
#       a. Summarize population
#       b. Find intersecting admin boundaries
#           i. For each admin boundary, clip the population raster, and summarize it
#######################################################################

import os, sys,inspect
import logging

import geopandas as gpd
import pandas as pd
import shapely
import rasterio
import fiona
from shapely.geometry import MultiPolygon
from rasterio.mask import mask
from rtree import index
from osgeo import gdal

cmd_folder = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)    
import GOSTRocks.misc
import GOSTRocks.rasterMisc

class compareAreas(object):
    def __init__(self, iso3, adminLayer, outputFolder='', popLayer='', urbanLayer='', loggingfile=''):
        if loggingfile == '':
            logging.basicConfig(level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
        else:
            logging.basicConfig(filename = loggingfile, level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
        
        logging.info("Starting urban admin comparison for %s" % iso3)
        
        self.urbanParams = GOSTRocks.misc.getUrbanParams()
        self.popLayer = self.urbanParams['lspopGrid'] if popLayer == '' else popLayer
        self.baseUrbanLayer = self.urbanParams['ls2012UrbanHDC_Vector'] if urbanLayer == '' else urbanLayer        
        self.outFolder = ("C:/Temp/%s_urban_inter/Try2" % iso3) if outputFolder == '' else outputFolder
        if not os.path.exists(self.outFolder):
            os.mkdir(self.outFolder)
        self.adminLayer = adminLayer
        self.iso3 = iso3 
        self.adminBoundaries = gpd.read_file(self.adminLayer)
        self.popStats = os.path.join(self.outFolder, "unioned_areas_population.csv")
        
        #Check CRS of input
        self.testCRS()
        
        #Generate Process Layers
        self.urbanLayer = self.getUrbanVector()        
    
    def testCRS(self):
        with rasterio.open(self.popLayer) as src:
            self.officialCRS = src.crs
                
        if self.adminBoundaries.crs != self.officialCRS:
            self.adminLayer = os.path.join(self.outFolder, "adminLayer.shp")
            if not os.path.exists(self.adminLayer):
                self.adminBoundaries = self.adminBoundaries.to_crs(self.officialCRS)
                self.adminBoundaries.to_file(driver = 'ESRI Shapefile', filename = self.adminLayer)
            
    def tabulateVectorResults(self, admCode, inputFile=''):
        ''' Combine the results from self.calculateVectorIntersection into simple urban classification
        
        methodology: http://ec.europa.eu/eurostat/web/rural-development/methodology
        
        [inputFile] - input unioned results, defaults to the results of self.calculateVectorIntersection
        RETURNS [ pandas df ] - tabulated results with EC classification
        '''
        if inputFile == '':
            unionedFile = self.popStats
        if not os.path.exists(unionedFile):
            loggine.warning("Input unioned file does not exist, run self.calculateVectorIntersection")
            raise ValueError("Missing input value")    
        unionRes = pd.read_csv(unionedFile)
        unionRes = unionRes.fillna(0)
        #Create UrbanPop Field
        unionRes['UrbanPop'] = unionRes['SUM'] * unionRes['GRIDCODE']
        unionRes['RuralPop'] = unionRes['SUM'] * abs(unionRes['GRIDCODE'] - 1)
        tabledRes = pd.pivot_table(unionRes, values=['UrbanPop','RuralPop'], index=[admCode])
        tabledRes['TotalPop'] = tabledRes['UrbanPop'] + tabledRes['RuralPop']
        tabledRes['PerRural'] = tabledRes['RuralPop'] / tabledRes['TotalPop']
        def define_urban(row):
            ''' http://ec.europa.eu/eurostat/web/rural-development/methodology
            '''
            if row['PerRural'] > 0.5:
                return "Rural"
            if row['PerRural'] < 0.5 and row['PerRural'] > 0.2:
                return "Intermediate"
            if row['PerRural'] < 0.2:
                return "Urban"

        tabledRes['UrbanClass'] = tabledRes.apply(lambda row: define_urban(row), axis=1)
        return tabledRes
        
    def getUrbanVector(self):
        '''Extract the urban areas from urban layer that intersect the admin areas
        '''
        urbanAdmin = os.path.join(self.outFolder, "%s_urbanAreas.shp" % self.iso3)
        if not os.path.exists(urbanAdmin):
            logging.info("Creating Urban Areas")
            inD = gpd.read_file(self.baseUrbanLayer)
            if inD.crs != self.officialCRS:
                logging.warning("CRS of urban data are being converted to CRS of admin file")
                inD = inD.to_crs(self.officialCRS)
            allAdmin = shapely.ops.cascaded_union(MultiPolygon(self.adminBoundaries['geometry'].tolist()))
            selUrban = GOSTRocks.misc.selectByIntersection(allAdmin, inD, exact=False)
            selUrban.crs = self.officialCRS
            selUrban.to_file(driver = 'ESRI Shapefile', filename = urbanAdmin)
            ###It would be interesting to see this as well, but not necessary right now
            #selUrban2 = GOSTRocks.misc.selectByIntersection(allAdmin, inD, exact=True)
            #selUrban2.to_file(driver = 'ESRI Shapefile', filename = urbanAdmin.replace(".shp", "_exact.shp"))
            logging.info("Created urban areas")
        return urbanAdmin
            
    def calculateVectorIntersection(self):
        ''' Union the input urban areas and input admin areas, then run zonal stats
            on the inputpopulation layer            
        '''
        #Create a union of the administrative boundaries and urban areas
        tempUnion = os.path.join(self.outFolder, "unioned_areas.shp")
        if not os.path.exists(tempUnion):
            logging.info("Creating union of urban areas and admin boundaries")
            urban = gpd.read_file(self.urbanLayer)
            if urban.crs != self.officialCRS:
                urban = urban.to_crs(self.officialCRS)
                
            urbanAdminUnion = gpd.overlay(urban, self.adminBoundaries, how='union')
            urbanAdminUnion.crs = self.officialCRS
            urbanAdminUnion.to_file(driver = 'ESRI Shapefile', filename = tempUnion)
        
        #Clip the population raster the the admin boundaries
        localPop = os.path.join(self.outFolder, "popLayer.tif")
        if not os.path.exists(localPop):     
            with rasterio.open(self.popLayer) as src:
                with fiona.open(self.adminLayer, "r") as shapefile:
                    geoms = [feature["geometry"] for feature in shapefile]                
                out_image, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta.copy()

            out_meta.update({"driver": "GTiff",
                             "height": out_image.shape[1],
                             "width": out_image.shape[2],
                             "transform": out_transform})

            with rasterio.open(localPop, "w", **out_meta) as dest:
                dest.write(out_image)

        #Run raster stats based on that union
        if not os.path.exists(self.popStats):
            logging.info("Running zonal statistics on unioned data")
            curRes = GOSTRocks.rasterMisc.zonalStats(tempUnion, localPop, bandNum=1, minVal = 0)
            zonalResults = pd.DataFrame(curRes, columns=['SUM','MIN','MAX','MEAN'])
            #bind the zonalstats to the boundaries
            adminBoundaries = gpd.read_file(tempUnion)
            zonalResults = pd.concat([adminBoundaries, zonalResults], axis=1)
            zonalShape = gpd.GeoDataFrame(zonalResults)
            zonalShape.to_file(driver="ESRI Shapefile", filename=self.popStats.replace(".csv", ".shp"))
            zonalShape.crs = self.officialCRS
            zonalResults = zonalResults.drop(['geometry'], axis=1)
            zonalResults.to_csv(self.popStats, header=True, index=True, encoding='utf-8')
        return self.popStats
