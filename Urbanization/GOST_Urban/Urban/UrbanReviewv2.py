###################################################################################################
# Urbanization Review
# Benjamin P. Stewart, Feb 2017
# Purpose: Calculate urbanization metrics for a defined country
#   1. Nighttime Lights
#		Extract footprints from DMSP data
#		Use Nighttime Lights urban extents for following summaries
#   2. GHSL
#		Map and summarize within NTL extents
#   3. GUF
#		Map and summarize within NTL extents
#   4. Gridded Population (Landscan 2012, as well as GHSPop
#		[Landscan 2012] Calculate urban extents following European Commission methodology
#		Summarize population with NTL extents
#	5. Calculate Puga Index on the GHSPop gridded population
#	6. Calculate Compactness and Discontiguity (from Harvard Toolbox) on GHSPop for all cities > 50000 people
#	7. Map NTL, GHSL, GUF, Population extents for every named NTL extent
#	8. TODO - Allow for extraction and zipping of data
#	9. TODO - Write all results to a single output excel file (for ease of data joining)
###################################################################################################
import os, sys, csv, arcpy, shutil, json, dbfpy, inspect
cmd_folder = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

#from dkan.client import DatasetAPI

from GOSTRocks.misc import *
from GOSTRocks.xlsxStuff import *
from GOSTRocks.arcpyMisc import *
from GOSTRocks.Urban.SummarizeGHSL import *
import GOSTRocks.metadataManagement
import GOSTRocks.Urban.NighttimeLightsFootprints
import GOSTRocks.Urban.NighttimeLightsThreshold
import GOSTRocks.Urban.WorldPop_Library as wp
import GOSTRocks.Urban.PugaIndex as puga
import GOSTRocks.Urban.MITUrbanForm.Compactness as MITCompactness
import GOSTRocks.Urban.MITUrbanForm.Discontiguity as MITDiscontiguity

def createFolder(bFolder, name):
    outF = os.path.join(bFolder, name)
    if not os.path.exists(outF):
        os.mkdir(outF)
    return (outF)

def calculateUrban(iso3, outputFolder, 
                   cNTL=True, cGHSL=True, cGUF=True, cPopulation=True, cAI=True, cAdmin=False, 
                   cPuga=True, cCompactness=True, cDiscontiguity=True, mapCities=True,
                   iUrban=False, finalizeProduct=False,
                   urbanExtents = "NTL", popIdx = "Pop", nameIdx = "Schnm",
                   popThresholds = [300,1500,5000,50000],
                   globalCities = "", admin0Polys = "", admin1Polys = "", admin2Polys = "",
                   ntlThresholds = "",urbanMetricsMetadata = "",gufMetadata = "",gufReferences = "",
                   gufTiles = "",gufFolder = "",gufMap = "",gufSymbology = "",lspopDensity = "",
                   lspopGrid = "",popMap = "",popHdStyle = "",popUrbStyle = "",landscanPopMetadata = "",
                   landscanDataMetadata = "",ghslFolder = "",ghsldatamaskFolder = "",ghslMetadata = "",
                   ghslOutline = "",ghslSymbology = "",dmLyrSymbology = "",extentsSymbology = "",
                   inNTLFolder = "",inNTLFolder_radCal = "",radCalT0 = "",radCalT1 = "",
                   ghspop75 = "",ghspop90 = "",ghspop00 = "",ghspop15 = "",ntlISA = "",aiUrban=""):
    '''
    
    '''
    #Read Parameters in from urban JSON parameters file
    jsonParams = getUrbanParams()
    
    globalCities = jsonParams['globalCities'] if globalCities == '' else globalCities 
    admin0Polys = jsonParams['admin0Polys'] if admin0Polys == '' else admin0Polys 
    admin1Polys = jsonParams['admin1Polys'] if admin1Polys == '' else admin1Polys
    admin2Polys = jsonParams['admin2Polys'] if admin2Polys == '' else admin2Polys
    ntlThresholds = jsonParams['ntlThresholds'] if ntlThresholds == '' else ntlThresholds
    urbanMetricsMetadata = jsonParams['urbanMetricsMetadata'] if urbanMetricsMetadata == '' else urbanMetricsMetadata
    gufMetadata = jsonParams['gufMetadata'] if gufMetadata == '' else gufMetadata
    gufReferences = jsonParams['gufReferences'] if gufReferences == '' else gufReferences
    gufTiles = jsonParams['gufTiles'] if gufTiles == '' else gufTiles
    gufFolder = jsonParams['gufFolder'] if gufFolder == '' else gufFolder
    gufMap = jsonParams['gufMap'] if gufMap == '' else gufMap
    gufSymbology = jsonParams['gufSymbology'] if gufSymbology == '' else gufSymbology
    lspopDensity = jsonParams['lspopDensity'] if lspopDensity == '' else lspopDensity
    lspopGrid = jsonParams['lspopGrid'] if lspopGrid == '' else lspopGrid
    popMap = jsonParams['popMap'] if popMap == '' else popMap
    popHdStyle = jsonParams['popHdStyle'] if popHdStyle == '' else popHdStyle
    popUrbStyle = jsonParams['popUrbStyle'] if popUrbStyle == '' else popUrbStyle
    landscanPopMetadata = jsonParams['landscanPopMetadata'] if landscanPopMetadata == '' else landscanPopMetadata
    landscanDataMetadata = jsonParams['landscanDataMetadata'] if landscanDataMetadata == '' else landscanDataMetadata
    ghslFolder = jsonParams['ghslFolder'] if ghslFolder == '' else ghslFolder
    ghsldatamaskFolder = jsonParams['ghsldatamaskFolder'] if ghsldatamaskFolder == '' else ghsldatamaskFolder
    ghslMetadata = jsonParams['ghslMetadata'] if ghslMetadata == '' else ghslMetadata
    ghslOutline = jsonParams['ghslOutline'] if ghslOutline == '' else ghslOutline
    ghslSymbology = jsonParams['ghslSymbology'] if ghslSymbology == '' else ghslSymbology
    dmLyrSymbology = jsonParams['dmLyrSymbology'] if dmLyrSymbology == '' else dmLyrSymbology
    extentsSymbology = jsonParams['extentsSymbology'] if extentsSymbology == '' else extentsSymbology
    inNTLFolder = jsonParams['NTLFolder'] if inNTLFolder == '' else inNTLFolder
    inNTLFolder_radCal = jsonParams['NTLFolder_radCal'] if inNTLFolder_radCal == '' else inNTLFolder_radCal
    radCalT0 = jsonParams['radCalT0'] if radCalT0 == '' else radCalT0
    radCalT1 = jsonParams['radCalT1'] if radCalT1 == '' else radCalT1
    ghspop75 = jsonParams['ghspop75'] if ghspop75 == '' else ghspop75
    ghspop90 = jsonParams['ghspop90'] if ghspop90 == '' else ghspop90
    ghspop00 = jsonParams['ghspop00'] if ghspop00 == '' else ghspop00
    ghspop15 = jsonParams['ghspop15'] if ghspop15 == '' else ghspop15
    ntlISA = jsonParams['ntlISA'] if ntlISA == '' else ntlISA    
    aiUrban = jsonParams['aiUrban'] if aiUrban == '' else aiUrban    
    
    #Create output folders and define output shapefiles
    baseOutputFolder = createFolder(outputFolder, "%s_v2" % iso3)  
    tempFolder = createFolder(baseOutputFolder, "Temp")    
    #2 Urban Geospatial Data
    urbanExtentFolder = createFolder(baseOutputFolder, "UrbanExtents")
    ntlFolder = createFolder(urbanExtentFolder, "NTL")
    popOutputFolder = createFolder(urbanExtentFolder, "EC_Pop")
    aiFolder = createFolder(urbanExtentFolder, "AI")    
    #3 Urban spatial indicators  Economic
    metricsFolder = createFolder(baseOutputFolder, "Economic_Metrics")
    #4 Urban spatial indicators  Morphological
    morphoFolder = createFolder(baseOutputFolder, "Morphological_Metrics")
    builtUpFolder = createFolder(morphoFolder, "BuiltUpAreas")
    analyticsFolder = createFolder(morphoFolder, "MIT_Metrics")
    #5Urban spatial indicators  Demographic
    populationOutputFolder = createFolder(baseOutputFolder, "Demographics_Metrics")
    #Other
    mapOutputFolder = createFolder(baseOutputFolder, "Maps")
    
    if finalizeProduct:
        zipOutputFolder = tempFolder
    
    #builtup area variables
    ghslShape = os.path.join(builtUpFolder, "%s_GHSL.shp" % iso3)
    GHSLoutputSummaryExcel = os.path.join(builtUpFolder, "%s_GHSL_Extent_Stats.xlsx" % iso3)    
    gufShape = os.path.join(builtUpFolder, "%s_GUF.shp" % iso3)
    tempGuf = os.path.join(tempFolder, "tempGUF.tif")
    aiUrbanAreas = os.path.join(aiFolder, "AI_Landscan2012.tif")
    builtAreaMetrics = os.path.join(builtUpFolder, "Morphological_Metrics.xlsx")
    
    #griddedPopulation variables
    tempPop = os.path.join(popOutputFolder, "populationGrid.tif")
    tempDen = os.path.join(popOutputFolder, "populationDen.tif")
    ghsPopStats = os.path.join(populationOutputFolder, "Population_Statistics.csv")
    urbClstGrid = os.path.join(popOutputFolder, "URB_CLST")
    hdClstGrid = os.path.join(popOutputFolder, "HD_CLST_RAW")
    hdClstGridSmooth = os.path.join(popOutputFolder, "hd_clst")
    hdClstGridSmoothBinary = os.path.join(popOutputFolder, "bhd_clst")
    urbClstTbl = os.path.join(popOutputFolder, "URB_UrbanizationSummary.dbf")
    hdClstTbl = os.path.join(popOutputFolder, "HD_UrbanizationSummary.dbf")
    urbShp = os.path.join(popOutputFolder, "URB_Extents.shp")
    hdShp = os.path.join(popOutputFolder, "HD_Extents.shp")
    
    #urban analytics variables
    compactnessShapefile = os.path.join(analyticsFolder, "MIT_compactness.shp")
    compactnessFile = os.path.join(analyticsFolder, "MIT_compactness.csv")
    discontiguityShapefile = os.path.join(analyticsFolder, "MIT_discontiguity.shp")
    urbanMetricsExcel = os.path.join(analyticsFolder, "%s_Urban_Metrics_Stats.xlsx" % iso3)
    pugaFile = os.path.join(analyticsFolder, "PUGA_Statistics.csv")
    
    #nighttime lights variables
    NTLoutputSummaryExcel = os.path.join(metricsFolder, "%s_NTL_Extent_Stats.xlsx" % iso3)
    finalShape = os.path.join(ntlFolder, "%s_masterExtents.shp" % iso3)
    mapShape = os.path.join(ntlFolder, "%s_masterExtents_forMapping.shp" % iso3)
    
    #Mapping variables
    outPopMap = os.path.join(mapOutputFolder, "population_Extent_Map.mxd")
    gufOutMap = os.path.join(mapOutputFolder, "GUF_Maps.mxd")   
    
    #This is the list of fields to not write out when writing shapefiles for urban metrics
    ignoreFields=["RC1996_T0","RC2010_T1","ntlChange","ntlChgCorr","intensive","extensive","extenCorr","areaChg","raw1992","raw1993","raw1994","raw1995","raw1996","raw1997","raw1998","raw1999","raw2000","raw2001","raw2002","raw2003","raw2004","raw2005","raw2006","raw2007","raw2008","raw2009","raw2010","raw2011","raw2012","rc1996","rc1999","rc2000","rc2002","rc2004","rc2005","rc2010","ExtTypeT1","CtyCntT1","ExtTypeT0","CtyCntT0"]
    
    #Create admin boundary files
    admin0 = os.path.join(tempFolder, "%s_0.shp" % iso3)
    admin1 = os.path.join(tempFolder, "%s_1.shp" % iso3)
    admin2 = os.path.join(tempFolder, "%s_2.shp" % iso3)
    if not arcpy.Exists(admin0):
        arcpy.Select_analysis(admin0Polys, admin0, '"ISO3" = \'%s\'' % iso3)
    if not arcpy.Exists(admin1):
        arcpy.Select_analysis(admin1Polys, admin1, '"ISO3" = \'%s\'' % iso3)
    if not arcpy.Exists(admin2):
        arcpy.Select_analysis(admin2Polys, admin2, '"ISO3" = \'%s\'' % iso3)

    #Calculate Nighttime lights footprints
    if cNTL:
        tPrint("***Processing NTL for %s" % iso3)
        #read the urban thresholds into a data dictionary
        countries = { rows[0]:rows[1:] for rows in csv.reader(open(ntlThresholds), delimiter=',') }
        newCities = os.path.join(ntlFolder, iso3 + os.path.basename(globalCities))    
        thresh = 0    
        try:
            thresh = countries[iso3][0]
        except: #If the threshold is not calculated, do so            
            thresh = GOSTRocks.Urban.NighttimeLightsThreshold.calculateThreshold(iso3, ntlThresholds)
        if not os.path.exists(finalShape):
            if not os.path.exists(ntlFolder):
                os.makedirs(ntlFolder)
        #Calculate NTL footprints
        if not arcpy.Exists(finalShape):
            GOSTRocks.Urban.NighttimeLightsFootprints.calculateFootprints(iso3, thresh, 
                ntlFolder, globalCities, radCalT0, radCalT1, inNTLFolder, inNTLFolder_radCal, admin0Polys,
                0.02, popIdx, nameIdx)        
        #Post Processing includes creating output tables, maps, and adding metadata to shapefiles
        GOSTRocks.Urban.NighttimeLightsFootprints.FPpostProcessing(metricsFolder,mapOutputFolder,ntlFolder,
            NTLoutputSummaryExcel, iso3, finalShape, newCities, thresh, popIdx, createMaps=False)            
    
    #Summarize GHSL within the defined urban footprints DEFAULT - Nighttime Lights Extents
    if cGHSL: 
        tPrint("***Processing GHSL for %s" % iso3)
        if urbanExtents == "NTL":
            arcpy.Select_analysis(finalShape, ghslShape, '"Pop" > 50000')
            #Delete most of the fields
            fieldsToDelete = []
            for f in arcpy.ListFields(ghslShape):
                if f.name in ignoreFields:
                    fieldsToDelete.append(f.name)
            arcpy.DeleteField_management(ghslShape, fieldsToDelete)
        else:            
            arcpy.CopyFeatures_management(urbanExtents, ghslShape)
        summarizeGHSL_toShp(ghslShape, ghslFolder, ghsldatamaskFolder, builtUpFolder, tempFolder, ghslOutline)
        writeShapefileXLS(builtAreaMetrics, ghslShape, "GHSL")
        
    #Summarize GUF within the defined urban footprints DEFAULT - Nighttime Lights Extents
    if cGUF:
        tPrint("***Processing GUF for %s" % iso3)
        if urbanExtents == "NTL":
            arcpy.Select_analysis(finalShape, gufShape, '"Pop" > 50000')
            #Delete most of the fields
            fieldsToDelete = []
            for f in arcpy.ListFields(gufShape):
                if f.name in ignoreFields:
                    fieldsToDelete.append(f.name)
            arcpy.DeleteField_management(gufShape, fieldsToDelete)
        else:            
            arcpy.CopyFeatures_management(urbanExtents, gufShape)
                        
        #Identify the intersecting GUF tiles
        summarizeGUF_toShp(gufShape, gufTiles, gufFolder, tempFolder)
        writeShapefileXLS(builtAreaMetrics, gufShape, "GUF")
        
    #Extract the Agglomeration Extents from the global dataset
    if cAI:
        arcpy.Clip_management(aiUrban, '#', aiUrbanAreas, admin0, clipping_geometry="ClippingGeometry")
        
    
    #Calculate the population based extents, based on EC classification methodology
    if cPopulation: 
        tPrint("***Processing gridded population for %s" % iso3)      
        lowDensVal = popThresholds[0]
        highDensVal = popThresholds[1]
        lowDensPop = popThresholds[2]
        highDensPop = popThresholds[3]
        
        #If the folder does not exist, create it 
        if not os.path.exists(popOutputFolder):
            os.mkdir(popOutputFolder)
        
        #Clip the population layers for processing
        if not arcpy.Exists(tempPop):
            arcpy.Clip_management(lspopGrid, "#", tempPop, admin0, '', "ClippingGeometry")
        if not arcpy.Exists(tempDen):
            arcpy.Clip_management(lspopDensity, "#", tempDen, admin0, '', "ClippingGeometry")
        
        #Create the high density and urban rasters if they do not exist
        tPrint("Create the high density and urban rasters if they do not exist")
        if not arcpy.Exists(hdClstGridSmooth):
            wp.createUrbanClusters(tempPop, tempDen, lowDensVal, urbClstGrid)
            wp.createUrbanClusters(tempPop, tempDen, highDensVal, hdClstGrid)
            tPrint("smoothing rasters")
            try:
                arcpy.env.workspace = popOutputFolder
                wp.smoothClusters(tempDen, tempPop, hdClstGrid, hdClstGridSmooth, highDensPop, 16, False)
                arcpy.CalculateStatistics_management(hdClstGridSmooth, "1", "1", "#", "OVERWRITE")
            except:
                tPrint("Threshold %s is not calculable as a high density cluster" % lowDensVal)
        
        #Summarize the output urban rasters - functions even if the rasters existed previously        
        arcpy.Delete_management(hdClstGrid)
        #Convert rasters to shapefiles
        if not arcpy.Exists(urbShp):
            extractFootprints(urbClstGrid, 4999, urbShp)
        if not arcpy.Exists(hdShp):
            extractFootprints(hdClstGridSmooth, 49999, hdShp)
        try:    
            straight = wp.summarizeWorldPop(urbClstGrid, tempPop, urbClstTbl, "ID")  
            smoothed = wp.summarizeWorldPop(hdClstGridSmooth, tempPop, hdClstTbl, "ID")              
        except:
            tPrint("Threshold %s does not have high density clusters" % lowDensVal)        
        
        #Run zonal statistics on GHSPop
        ghsPopRasters = [ghspop75, ghspop90, ghspop00, ghspop15, lspopGrid]
        if urbanExtents == "NTL":
            popExtents = finalShape			
        else:            
            popExtents = urbanExtents        
        popRes = summarizeGlobalData(popExtents, "FID", ghsPopRasters, type = ['N'] * len(ghsPopRasters))
        writeDict(popRes["Results"], ghsPopStats, popRes["Titles"])   
        
        
        
    if cPuga:
        tPrint("*****Calculating Puga Index")
        #Calculate the puga index
        if urbanExtents == "NTL":
            pugaIndex = finalShape
        else:            
            pugaIndex = urbanExtents
        if not arcpy.Exists(pugaFile):
            puga.main(pugaIndex, pugaFile, tempFolder, ghspop75, ghspop90, ghspop00, ghspop15)
        
        #Write results to output shapefiles
        writeShapefileXLS(urbanMetricsExcel, pugaIndex, "pugaExtents", ignoreFields=ignoreFields)
        writeCSVXLS(urbanMetricsExcel, pugaFile, "pugaValues")
            
            
    if cCompactness:
        tPrint("*****Calculating Compactness")
        #Compactess analyses individual points or polygons within a single city - in this case, we want to run it individually
        #   On every city larger than 50000 people on each of 4 ghsPop layers
        if not arcpy.Exists(compactnessShapefile):
            arcpy.Select_analysis(finalShape, compactnessShapefile, '"Pop" > 50000')
            compactness75 = MITCompactness.processUrbanExtents(compactnessShapefile, ghspop75, tempFolder, "p75")
            tPrint("Calculated Compactness for 1975")
            compactness90 = MITCompactness.processUrbanExtents(compactnessShapefile, ghspop90, tempFolder, "p90")
            tPrint("Calculated Compactness for 1990")
            compactness00 = MITCompactness.processUrbanExtents(compactnessShapefile, ghspop00, tempFolder, "p00")
            tPrint("Calculated Compactness for 2000")
            compactness15 = MITCompactness.processUrbanExtents(compactnessShapefile, ghspop15, tempFolder, "p15")
            tPrint("Calculated Compactness for 2015")

        #Write results to output excel table
        writeShapefileXLS(urbanMetricsExcel, compactnessShapefile, "Compactness", sortField="pop", sortOrder="D", ignoreFields=ignoreFields)
    
    if cDiscontiguity:
        tPrint("*****Calculating MIT Discontiguity")
        if not arcpy.Exists(discontiguityShapefile):
            arcpy.Select_analysis(finalShape, discontiguityShapefile, '"Pop" > 50000')
            curTiles = selectGHSLtiles(discontiguityShapefile, ghslOutline, 'Path', False) 
            if not arcpy.Exists(os.path.join(tempFolder, "ghslMaster.tif")):
                createSingleLayer(curTiles, tempFolder, "ghslMaster.tif")
            curRasters = createYearlyLayers(os.path.join(tempFolder, "ghslMaster.tif"), tempFolder)
            gufIdx = 0
            print curRasters
            for cTile in curRasters:
                gufIdx = gufIdx + 1
                MITDiscontiguity.calcDiscontiguityFromGUF(discontiguityShapefile, cTile, tempFolder, "GHSL%s" % gufIdx)
            
        #Write results to output excel table       
        writeShapefileXLS(urbanMetricsExcel, discontiguityShapefile, "Discontiguity", sortField="pop", sortOrder="D", ignoreFields = ignoreFields)
        
    if mapCities:
        tPrint("Creating city maps")
        if urbanExtents == "NTL":
            arcpy.Select_analysis(finalShape, mapShape, '"Pop" > 50000')
            #Delete most of the fields
            fieldsToDelete = []
            for f in arcpy.ListFields(mapShape):
                if not f.name in ["ExtentName","gAreaKM","pop","FID","Shape"]:
                    fieldsToDelete.append(f.name)
            arcpy.DeleteField_management(mapShape, fieldsToDelete)
        else:            
            arcpy.CopyFeatures_management(urbanExtents, mapShape)

        #Population maps
        wp.mapPop(popMap, outPopMap, hdClstGridSmooth, urbClstGrid, popHdStyle, popUrbStyle, mapShape, extentsSymbology)
        createMapFromMXD_loop(arcpy.mapping.MapDocument(outPopMap), os.path.join(mapOutputFolder, "Population_Map.png"), "%s_masterExtents_forMapping" % iso3, "ExtentName")
        
        #GUF Maps
        gufTiles = selectGHSLtiles(gufShape, gufTiles, 'FileName', False)        
        addGUF(gufMap, gufTiles, gufShape, gufOutMap, gufFolder, gufSymbology, extentsSymbology)
        createMapFromMXD_loop(arcpy.mapping.MapDocument(gufOutMap), os.path.join(mapOutputFolder, "GUF_Map.png"), "%s_GUF" % iso3, "ExtentName")

        #Map GHSL       
        ghslMap = os.path.join(mapOutputFolder, "GHSL_Maps.mxd")
        addGHSL(mapShape, ghslMap, ghslFolder, ghsldatamaskFolder, ghslSymbology, extentsSymbology, ghslOutline)
        createMapFromMXD(arcpy.mapping.MapDocument(ghslMap), os.path.join(mapOutputFolder, "GHSL_Map.png"))
        createMapFromMXD_loop(arcpy.mapping.MapDocument(ghslMap), os.path.join(mapOutputFolder, "GHSL_Map.png"), "%s_masterExtents_forMapping" % iso3, "ExtentName")

        #Nighttime Lights
        createMapFromMXD_loop(arcpy.mapping.MapDocument(os.path.join(mapOutputFolder, "%s NighttimeLights_Extents.mxd" % iso3)),
            os.path.join(mapOutputFolder, "NTL_Map.png"), "Master Footprints", "ExtentName",
            statusDefinition = {'Master Footprints':"NOT \"ExtentName\" = '0'"})        
            
    if iUrban:            
        tPrint("Calculating iUrban")
        #Variables for iUrban work
        iUrbanConfig = os.path.join(tempFolder, "%s_iUrban.json" % iso3)
        tempPop = os.path.join(popOutputFolder, "populationGrid.tif")
        tempDen = os.path.join(popOutputFolder, "populationDen.tif")
        cntryISA =  os.path.join(popOutputFolder, "ntlISA_%s.tif" % iso3)
        cntryGUF =  os.path.join(popOutputFolder, "guf_%s.tif" % iso3)
        cntryRoads =  os.path.join(popOutputFolder, "roads_%s.tif" % iso3)    
        #Generate the iUrban JSON configuration file
        iUrbanDict = {"country": iso3,
            "population_raster": tempPop,
            "population_dataset_type": "LandScan",
            "reference_year": 2012,
            "main_urban_patch_method": {
                "type": "default",
                "urban_patch_lat": 19.790259,
                "urban_patch_lon": 96.130120,
                "metro_threshold": 71712,
                "comment": "Capital"
            },
            "isa_dataset": cntryISA,
            "bua_dataset": tempGuf,
            "grip_data": cntryRoads,
            "administrative_level0": admin0,
            "administrative_level1": admin1,
            "inventory_pager_file": "data/asset_value_calc_files/cri_inventory_pager.csv",
            "inventory_popweight_file": "data/asset_value_calc_files/cri_inventory_popweight.csv",
            "ppd_prov_file": "data/asset_value_calc_files/cri_ppd_prov.csv",
            "dwelling_size_file": "data/asset_value_calc_files/cri_dwelling_size.csv",
            "res_ucc_file": "data/asset_value_calc_files/cri_res_ucc.csv",
            "nonres_exp_file": "data/asset_value_calc_files/cri_nonres_exp.csv",
            "nonres_occupancy_file": "data/asset_value_calc_files/cri_inventory_nonres.csv",
            "nonres_ucc_file": "data/asset_value_calc_files/cri_nonres_ucc.csv",
            "comments": "urbanratio (2010 UN-WUP), urbanratio_ref (2012 UN-WUP), pop (2010 UN WDI), pop_ref (2012 UN WDI)"
        }
        #Write the above json to the json file
        with open(iUrbanConfig, 'w') as outFile:
            json.dump(iUrbanDict, outFile)
        
        #Merge the guf tiles into one image
        if not arcpy.Exists(tempGuf):
            gufTiles = selectGHSLtiles(admin0, gufTiles, 'FileName', False)        
            fullPathGuf = []
            for g in gufTiles:
                fullPathGuf.append(os.path.join(gufFolder, g))      
            print ";".join(fullPathGuf)
            arcpy.MosaicToNewRaster_management(";".join(fullPathGuf), tempFolder, "tempGuf.tif", number_of_bands=1)        
        
        #runUrban.main([iUrbanConfig, tempFolder, "population"])        
        command = '"%s/wb_cdrp/model/run_all.py" %s %s --process_until population' % (cmd_folder, iUrbanConfig, tempFolder)
        print command
        os.system(command)                
        
    if finalizeProduct:
        #Upload the dataset to the New Data Catalog
        #Input Variables
        ddhurl = jsonParams['ddhurl']
        user = jsonParams['ddhUser']
        password = jsonParams['ddhPassword']
        api = DatasetAPI(ddhurl, user, password, False)
        #Process the list values dataset
        taxonomyLink = "%s%s" % (ddhurl, "/internal/listvalues")
        taxResults = api.get(taxonomyLink)
        countryResults = {}
        countryTypes = {}
        for entry in taxResults.json():
            try:
                if entry['machine_name'] == "field_wbddh_country":
                    countryResults[entry["list_value_name"]] = entry["tid"]        
                if entry['machine_name'] == "field_wbddh_economy_coverage":
                    countryTypes[entry["list_value_name"]] = entry["tid"]
            except:
                FUBAR = "Do Fucking Nothing"
        #Read in the admin0 shapefile to get the country code and the lending category codes
        taxonomyLink = "%s%s" % (ddhurl, "/internal/listvalues")
        taxResults = api.get(taxonomyLink)
        countryResults = {}
        countryTypes = {}
        for entry in taxResults.json():
            try:
                if entry['machine_name'] == "field_wbddh_country":
                    countryResults[entry["list_value_name"]] = entry["tid"]        
                if entry['machine_name'] == "field_wbddh_economy_coverage":
                    countryTypes[entry["list_value_name"]] = entry["tid"]
            except:
                FUBAR = "Do Fucking Nothing"
        db = dbf.Dbf(admin0.replace(".shp",".dbf"))
                
        countryCode = countryResults[db[0]["WB_ADM0_NA"]]#['Myanmar']#
        countryCategory = countryTypes[db[0]["LENDINGC"]]#['IDA']#
        datasetData = {
            'type':'dataset',
            "field_contact_email":"gost@worldbank.org",    
            "field_wbddh_data_class": {"und":{"id":"358"}},             #Data are public
            "field_wbddh_data_type": {"und":{"id":"295"}},              #Geospatial Data
            "field_frequency": {"und":"18"},
            "field_wbddh_terms_of_use": {"und":{"id":"434"}},           #Open Data Access
            "field_wbddh_languages_supported": {"und":{"id":"337"}},    #English    
            "field_granularity_list": {"und":{"id":"939"}},
            "field_topic":{"und":{"id":"384"}},
            "field_wbddh_working_unit_user": {"und":{"id":"1003"}},     #???
            "workflow_status":"published",
            "body":{"und":[{"value":"Urbanization Metrics for %s, describing urban shape, extent, and other characteristics. \n\n \
                Teams across the World Bank have worked with indices and tools which use urban spatial data to analyze and \
                characterize various aspects of urbanization, to facilitate the formulation of policies or the design of projects. \
                These often take advantage of the increasing profusion of spatial data sources, particularly the recent availability \
                of global data which allow comparison across cities and countries. \n\n The lack of shared information about \
                these tools has sometimes led to duplication of effort and lack of comparable outputs. To remedy this, \
                this data package attempts to summarize the tools currently in use for the reference of teams seeking to \
                carry out some form of urban spatial analysis. It builds on an accompanying document which discusses sources of urban spatial data." % db[0]["WB_ADM0_NA"]}]},    
            
            'title':'%s Urban Metrics' % db[0]["WB_ADM0_NA"],
            "field_wbddh_country": {"und":{"id":"%s" % countryCode}},                 #Afghanistan
            "field_wbddh_economy_coverage": {"und":{"id":"%s" % countryCategory}},       #WBG Economy Classification - Blend vs IDA
            
            "field_wbddh_dsttl_upi": {"und":{"autocomplete_hidden_value":"45981"}},      #ID of TTL - this is Ben
            "field_wbddh_map_projection":"",
            "field_wbddh_limitations_geo":"Benjamin Stewart/GOST is not necessarily the originator of this data, if you have questions about this data please contact gost@worldbank.org"           
        }
        dataset = api.node('create', data=datasetData)
        dataset_nid = dataset.json()['nid']
        print ("%s has been uploaded as %s" % (iso3, dataset_nid))
        
        #Final steps include deleting the temp folder, zipping all the results folders to DDH drive, and creating the DDH datasets
        resourceList = []
        resources = [populationOutputFolder,metricsFolder,morphoFolder,mapOutputFolder,urbanExtentFolder]
        resourceTitles = ["Demographics", "Economic Metrics", "Morphological Metrics", "Urban Maps", "Geospatial Data on Urban Extents"]
        resourceDescriptions = ["Population metrics describe the population derived from GHSPop and Landscan 2012 within the urban extents calculated by nighttime lights",
            "Economic metrics describe the change in nighttime lights for urban areas, which are a common proxy for assessing economic development",
            "Morphological metrics focus on the change in built-area of the urban extents, and calculate a number of urban metrics derived from an existing toolbox developed at MIT, along with another metric called Puga",
            "The maps folder contains a series of maps from all the various metrics for the largest cities in %s" % db[0]["WB_ADM0_NA"],
            "The urban extents define the urban areas across the country measured three ways: nighttime lights, the european cluster methodology run against Landscan 2012, and the agglomeration index"]
        for idx in range(0,5):  #create zip files 
            zipFolder = resources[idx]
            outZipFile = os.path.join(zipOutputFolder, os.path.basename(os.path.normpath(zipFolder)))
            if not os.path.exists("%s.zip" % outZipFile):
                xx = shutil.make_archive(outZipFile, 'zip', zipFolder)
            resourceData = {
                'type':'resource',
                'title':resourceTitles[idx], 
                "field_wbddh_data_class": {"und":{"id":"358"}},             #Data are public   
                "workflow_status":"published",                          
                'field_dataset_ref': {'und':{'target_id': dataset_nid}},
                'field_wbddh_resource_type' : {'und':{'id':'443'}},         #Resource is a dataset; documentation = 985
                'field_format' : {'und':{'id':'957'}},                       #Resource is a shp zip;
                
                "body":{"und":[{"value":resourceDescriptions[idx]}]}
            }

            resource = api.node('create', data=resourceData)
            resource_nid = resource.json()['nid']
            r = api.attach_file_to_node("%s.zip" % outZipFile, resource_nid, 'field_upload')
            resourceUpdate = "%s (%s)" % (resourceTitle, resource_nid)
            resourceList.append({'target_id':resourceUpdate})
            datasetWithResource = api.node('update', node_id=dataset_nid, data={"field_resources":{"und":resourceList}})
            print "%s has been uploaded as %s" % (os.path.basename(os.path.normpath(zipFolder)), resource_nid)
        
        shutil.rmtree(tempFolder)   #delete temp folder        
