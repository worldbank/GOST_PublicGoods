###############################################################################
# Extract Nighttime Lights footprints
# Benjamin P Stewart, April 2015
# Purpose: Functions to extract nighttime lights urban footprints
###############################################################################
import arcpy, os, datetime, csv, sys, re, openpyxl, shutil, tempfile #xlsxwriter
from dbfpy import dbf
from arcpy.sa import *
sys.path.append(r"C:\Users\wb411133\Box Sync\AAA_BPS\Code\GOST")
from GOSTRocks.arcpyMisc import *
from GOSTRocks.misc import *
from GOSTRocks.xlsxStuff import *
import GOSTRocks.metadataManagement

tempFolder = r"C:\Temp"
#Fields that define city distance fields
t0FID = "nearT0FID"
t1FID = "nearT1FID" 
t0Dist = "t0Dist"
t1Dist = "t1Dist"
strFields = ["ExtentType", "ExtentName", "Status"]
numFields = ["cityCount","pop","RC1996_T0","RC2010_T1","ntlChange","ntlChgCorr","intensive","extensive","extenCorr","areaChg",
    "raw1992","raw1993","raw1994","raw1995","raw1996","raw1997","raw1998","raw1999","raw2000",
    "raw2001","raw2002","raw2003","raw2004","raw2005","raw2006","raw2007","raw2008","raw2009",
    "raw2010","raw2011","raw2012","rc1996","rc1999","rc2000","rc2002","rc2004","rc2005","rc2010"]
fields = strFields + numFields
fpIdField = "FID"
deleteFiles = []
joinFields = ["CityType", "AggName"]
cLyr = 'cities'


###Convert Nighttime lights data into vector urban footprints
#countryLyr - input country from which we are extracting footprints 
    #####BEN: This may cause problems with cross border footprints, eh#####
#NTLFile - input NTL raster data
#threshold - cutoff for NTLFile. Above is urban, below is not
#outPolygons - output vector data file
def extractFootprints(countryLyr, NTLFile, threshold, outPolygons, tempFolder = r"C:\Temp"):
    tempRaster = os.path.join(tempFolder, "cntryRaster.tif")
    tempShp = os.path.join(tempFolder, "cntry.shp")
    tempShp2 = os.path.join(tempFolder, "cntry2.shp")
    tempShp3 = os.path.join(tempFolder, "cntry3.shp")
                                    
    #extract NTL values for current country
    if not countryLyr == -1:
        arcpy.Clip_management(NTLFile, '#', tempRaster, countryLyr, '#', 'ClippingGeometry', '#')
        tRaster = Raster(tempRaster)
    else:
        tRaster = Raster(NTLFile)
    tRaster = tRaster > float(threshold)

    #Convert to vector and select out the non zero values
    arcpy.RasterToPolygon_conversion(tRaster, tempShp, 'NO_SIMPLIFY', '#')
    arcpy.Select_analysis(tempShp, tempShp2, ' "GRIDCODE" > 0 ')
    arcpy.Buffer_analysis(tempShp2, tempShp3, "5 Meters")#, "FULL", "ROUND", "LIST", "GRIDCODE")
    arcpy.Dissolve_management(in_features=tempShp3, out_feature_class=outPolygons, dissolve_field="GRIDCODE", statistics_fields="", multi_part="SINGLE_PART", unsplit_lines="DISSOLVE_LINES")
    #Add area field and calculate area in square kilometres
    tryAddField(outPolygons, "gAreaKM", "FLOAT")
    arcpy.CalculateField_management(outPolygons, "gAreaKM", "!SHAPE.AREA@SQUAREKILOMETERS!", "PYTHON")
    for f in [tempRaster, tempShp, tempShp2, tempShp3]:
        arcpy.Delete_management(f)

# Run zonal statistics on all nighttime lights data for inMaskData    
def summarizeAllNTL(inMaskData, idField, oFolder, version, inNTLFolder, inNTLFolder_radCal):
    outFiles = []
    for yr in range(1992, 2013, 1):
        #Get raster vairable - This step will vary based on how your structure is set-up
        #This current version works if all the images are in a single folder defined in inNTLFolder
        baseFolder = os.path.join(inNTLFolder, str(yr))
        ntlFiles = os.listdir(baseFolder)        
        for f in ntlFiles:
            if re.search("ElvidgeCorrected.tif$", f):
                ntlFile = os.path.join(baseFolder, f)
                outTable = os.path.join(oFolder, "%s_%s_%s.dbf" %("NTLZonal",yr,version))
                outFiles.append(outTable)
                ZonalStatisticsAsTable (inMaskData, idField, ntlFile, outTable)
    #Run again on all radiance calibrated data
    inFiles = sorted(os.listdir(inNTLFolder_radCal)) 
    for f in inFiles:
        if re.search("avg_vis_Corrected.tif$", f):
            ntlFile = os.path.join(inNTLFolder_radCal, f)
            outTable = os.path.join(oFolder, "%s_%s_%s.dbf" %("NTLZonal",f.replace(".",""),version))
            outFiles.append(outTable)
            ZonalStatisticsAsTable (inMaskData, idField, ntlFile, outTable)
    allZonal = summarizeDbf(outFiles, "FID_")
    deleteFiles.extend(outFiles)    
    inMaskData = None
    return allZonal

def calculateNearFP(inShpT0, inShpT1, newCities):  
    #Add fields to cities shapefile for use in calculating nearest footprints
    inFields = { f.name for f in arcpy.ListFields(newCities) }     
    if not t0FID in inFields:
        tryAddField(newCities, t0FID, "LONG")
        tryAddField(newCities, t1FID, "LONG")
        tryAddField(newCities, t0Dist, "FLOAT")
        tryAddField(newCities, t1Dist, "FLOAT")
    
    #For each feature in the input cities, calculate distance to the nearest footprint and rename the field   
    arcpy.Near_analysis(newCities, inShpT0)
    arcpy.CalculateField_management (newCities, t0FID, "!NEAR_FID!", "PYTHON_9.3") 
    arcpy.CalculateField_management (newCities, t0Dist, "!NEAR_DIST!", "PYTHON_9.3") 
    arcpy.Near_analysis(newCities, inShpT1)
    arcpy.CalculateField_management (newCities, t1FID, "!NEAR_FID!", "PYTHON_9.3") 
    arcpy.CalculateField_management (newCities, t1Dist, "!NEAR_DIST!", "PYTHON_9.3") 
    inShpT0 = None
    inShpT1 = None
    
#Get the summed population
def getCityPopulation(cityLyr, cityCnt, searchDist, popIdx, nameIdx):
    cityCursor = arcpy.SearchCursor(cityLyr)
    pop = 0
    maxPop = 0
    outCity = ''
    fieldNames = { f.name for f in arcpy.ListFields(cityLyr) }
    cty = cityCursor.next()
    cFid = 0
    while cty:
        if cityCnt == 1:
            #When searching for the nearest t0footprint, need to only select features 
            #   that are within the defined search distance            
            cty0Dist = cty.getValue(t0Dist)
            if cty0Dist > searchDist:
                cty0ID = -1
            else:
                cty0ID = cty.getValue(t0FID)
            
            return([cty.getValue("FID"), cityCnt, 'Single City', cty.getValue(nameIdx), 
                cty.getValue(popIdx), cty0ID])
        else:
            ctyPop = cty.getValue(popIdx)
            pop = pop + float(ctyPop)
            t0FIDs = (cty.getValue(t0FID))
            if ctyPop > maxPop:
                maxPop = ctyPop
                #Get the city name and FID for storing in the agglomeration database
                outCity = cty.getValue(nameIdx)
                cFid = cty.getValue("FID")
        cty = cityCursor.next()
    if cityCnt > 1:
        return([cFid, cityCnt, 'Agglomeration', outCity, pop, t0FIDs])   

def createMasterFP(t0, t1, finalShp, searchDist):
    #Update t1 with values from t0 - ntlChange values and ntlT0 type values
    tryAddField(t1, "ExtTypeT0", "TEXT")
    tryAddField(t1, "CtyCntT0")
    arcpy.MakeFeatureLayer_management(t0, "t0")
    arcpy.CalculateField_management(t1, "Status", "APPEAR")
    fpCursor = arcpy.da.UpdateCursor(t1, ["FID","Status","ExtTypeT0", "CtyCntT0","RC1996_T0","RC2010_T1","ntlChange","ntlChgCorr","intensive","extensive","extenCorr","areaChg","rc1996","rc2010","gAreaKM","ExtentName"])
    for row in fpCursor:
        #Search for features in T0 that intersect the current shape
        arcpy.MakeFeatureLayer_management(t1, "t1", '"FID" =  %s' % row[0])
        arcpy.SelectLayerByLocation_management("t0", "WITHIN_A_DISTANCE", "t1", 0)
        t0Cursor = arcpy.da.SearchCursor("t0", ["ExtTypeT0", "CtyCntT0", "rc1996", "rc2010","gAreaKM","ExtentName"])
        RCT01996 = 0 
        RCT02010 = 0
        RCT11996 = row[12]
        RCT12010 = row[13]
        t0Count = 0
        tArea = 0
        t0Type = ""
        status = "APPEAR"
        for rowT0 in t0Cursor:
            status = "FOUND"
            RCT01996 += rowT0[2]
            RCT02010 += rowT0[3]
            tArea += rowT0[4]        
            if t0Type != "Agglomeration":
                t0Type = rowT0[0]
                if t0Type != "Single City":                
                    t0Type = rowT0[0]
            if rowT0[1] > 0:
                t0Count += rowT0[1]
        ntlChangeCorrected = RCT12010-RCT01996 - (RCT11996 - RCT01996)
        extensiveCorrected = (RCT12010 - RCT02010) - (RCT11996 - RCT01996)
        #          ["FID",  "status" "CityTypeT0",  "CityCountT0","t0NTL",  "t1NTL"
        row[0:6] = [row[0], status,   t0Type,        t0Count,      RCT01996, RCT12010]
        #           "ntlChange",        "ntlChangeCorrected", "intensive",          "extensive",          "extensiveCorrected", "areaChg"
        row[6:12] = [RCT12010-RCT01996,  ntlChangeCorrected,   RCT02010 - RCT01996, (RCT12010 - RCT02010), extensiveCorrected,  row[14] - tArea]
        fpCursor.updateRow(row)

    #Need to include DISAPPEARING features from T0 in the final Master Footrint file
    arcpy.MakeFeatureLayer_management(t1, "t1All")
    arcpy.SelectLayerByAttribute_management('t0', 'CLEAR_SELECTION')
    arcpy.SelectLayerByAttribute_management('t1All', 'CLEAR_SELECTION')
    arcpy.SelectLayerByLocation_management('t0', "INTERSECT", 't1All')
    arcpy.SelectLayerByAttribute_management ("t0", "SWITCH_SELECTION")
    arcpy.CalculateField_management("t0", "Status", '"DISAPPEAR"')
    arcpy.Merge_management(["t1All","t0"], finalShp)
    
def attributeT0NTL(t0NTL, inFP, searchDist, popIdx, nameIdx, timePeriod):
    #Add all values to the in footprints
    upCursor = arcpy.UpdateCursor(inFP)
    row = upCursor.next()
    fidIdx = t1FID
    distIdx = t1Dist
    missedVals = "MISSED"
    #Search for different values in the cities file depending on the time period
    if timePeriod == 0:
        fidIdx = t0FID
        distIdx = t0Dist
        missedVals = "MISSED"
    while row:
        allVals = [0] * 13    
        query = '%s = %s AND %s < %s' % (fidIdx, row.getValue("FID"), distIdx, searchDist)            
        arcpy.SelectLayerByAttribute_management("cityLyr", "NEW_SELECTION", query)
        t0CityIntersect = int(arcpy.GetCount_management('cityLyr').getOutput(0))               
        if t0CityIntersect > 0:
            t1CityInfo = getCityPopulation('cityLyr', t0CityIntersect, searchDist, popIdx, nameIdx)                         
            allVals[0:5] = [t1CityInfo[2],t1CityInfo[3], "FOUND", t1CityInfo[1], t1CityInfo[4]]
        else:
            allVals[2] = missedVals           
            
        allVals.extend(t0NTL[row.getValue("FID")])
        for attrIdx in range(0, len(allVals)):
            row.setValue(fields[attrIdx], allVals[attrIdx])            
        upCursor.updateRow(row)
        row = upCursor.next()

def addFields(inShp):
    fieldNames = { f.name for f in arcpy.ListFields(inShp) }
    #Add all string and numeric fields to the feature class (These should be doable on one line, eh?)
    for fName in strFields:
        if not fName in fieldNames:
            tryAddField(inShp, fName, "TEXT") 
    for fName in numFields:
        if not fName in fieldNames:
            tryAddField(inShp, fName, "FLOAT")  
    inShp = None

def calculateFootprints(cntry, thresh, outputFolder, inputCities, 
        radCalT0, radCalT1, inNTLFolder, inNTLFolder_radCal, inputCountries,
        searchDist=0.02, popIdx = "ES00POP", nameIdx = "SCHNM"):
    
    newCities = os.path.join(outputFolder, cntry + os.path.basename(inputCities))
    finalShape = os.path.join(outputFolder, "%s_masterExtents.shp" % cntry)
    fpT0 = os.path.join(outputFolder, "%s_radCalT0.shp" % cntry)
    fpT1 = os.path.join(outputFolder, "%s_radCalT1.shp" % cntry)
    #Create a feature layer for the current country
    cntryLyr = arcpy.MakeFeatureLayer_management (inputCountries,"temp_lyr", "ISO3 = '%s'" % str(cntry))
    #Create a new city layer for just the active country
    arcpy.Clip_analysis(inputCities, cntryLyr, newCities)    
    arcpy.MakeFeatureLayer_management(newCities, 'cityLyr')          
    if arcpy.Exists(fpT0) and arcpy.Exists(fpT1):
        tPrint("%s already exists" % (cntry))
    else:
        tPrint("Extracting footprints for %s" % cntry)    
        extractFootprints(cntryLyr, radCalT0, thresh, fpT0)
        tPrint("Processed RadCal T0: %s" % thresh)
        extractFootprints(cntryLyr, radCalT1, thresh, fpT1)
        tPrint("Processed RadCal T1")

    #Add fields to the T1 footprints file
    addFields(fpT1)
    addFields(fpT0)
    #Calculate distance from cities to footprints for joining cities to footprints
    calculateNearFP(fpT0, fpT1, newCities)
    
    #Run zonal statistics on the newly minted nighttime lights urban footprints
    #Specific subsets of the overall zonal statistics are extracted for use in the next procedures
    zonalNTL = summarizeAllNTL(fpT1, fpIdField, outputFolder,"T1", inNTLFolder, inNTLFolder_radCal)    
    fpZonalT1 = {}
    for i in range(0, len(zonalNTL)):
        curZonal = []        
        for j in [0,0,0,0,-1,-1,-1,-1]:            
            curZonal.append(zonalNTL[i][j])
        fpZonalT1[i] = curZonal

    zonalNTLT0 = summarizeAllNTL(fpT0, fpIdField, outputFolder,"T0", inNTLFolder, inNTLFolder_radCal)    
    fpZonalT0 = {}
    for i in range(0, len(zonalNTLT0)):
        curZonal = []
        for j in [0,0,0,0,-1,-1,-1,-1]:            
            curZonal.append(zonalNTLT0[i][j])
        fpZonalT0[i] = curZonal
        
    tPrint("Ran All nighttime lights zonal statistics on %s" % cntry)

    #Merge footprint values with city values and calculate night light change values    
    attributeT0NTL(zonalNTLT0, fpT0, searchDist, popIdx, nameIdx, 0)
    renameField(fpT0, "ExtentType", "ExtTypeT0")
    renameField(fpT0, "cityCount", "CtyCntT0")   
    attributeT0NTL(zonalNTL, fpT1, searchDist, popIdx, nameIdx, 1)
    renameField(fpT1, "ExtentType", "ExtTypeT1")
    renameField(fpT1, "cityCount", "CtyCntT1")   
    
    #Calculate change values in a master footprint file
    createMasterFP(fpT0, fpT1, finalShape, searchDist)
    tPrint("Merged City footprints with city data")
    
    #Join the nighttime lights footprints names and agglomerations definitions to the city shapefile        
    arcpy.JoinField_management(newCities, t1FID, fpT1, "FID", ["ExtTypeT1", "CtyCntT1", "Status", "ExtentName"])
    arcpy.MakeFeatureLayer_management(newCities, cLyr)
    arcpy.SelectLayerByAttribute_management(cLyr, "NEW_SELECTION", '"%s" > %s' %(t1Dist,searchDist))
    arcpy.CalculateField_management(cLyr, "ExtTypeT1", '""')      
    arcpy.CalculateField_management(cLyr, "Status", '""')     
    arcpy.CalculateField_management(cLyr, "ExtentName", '""')     
    arcpy.CalculateField_management(cLyr, "CtyCntT1", '0')      
    
    #Join the T0 information, calculate blanks for the empty features
    arcpy.JoinField_management(newCities, t0FID, fpT0, "FID", ["ExtTypeT0", "CtyCntT0", 'ExtentName'])   
    renameField(newCities, "ExtentNa_1", "ExtentT0")
    arcpy.SelectLayerByAttribute_management(cLyr, "NEW_SELECTION", '"%s" > %s' %(t0Dist,searchDist))
    arcpy.CalculateField_management(cLyr, "ExtTypeT0", '""')      
    arcpy.CalculateField_management(cLyr, "CtyCntT0", '0')      
    arcpy.CalculateField_management(cLyr, "ExtentT0", '""')      
    
    arcpy.SelectLayerByAttribute_management(cLyr, "NEW_SELECTION", '"CtyCntT0" > 0 AND "CtyCntT1" = 0')
    arcpy.CalculateField_management(cLyr, "Status", '"DISAPPEAR"')      
    arcpy.CalculateField_management(cLyr, "ExtentName", '""')      
            
    #Cleanup
    cntryLyr = None    
    for f in deleteFiles:
        arcpy.Delete_management(f)  

#For the output nighttime lights mxd fix the reference to the nighttime lights footprints
def fixMapConnections(mxd, df, cntry, outFolder):
    for lyrName in ['Rad Cal Footprint 2010', 'Rad Cal Footprint 1996', 'Master Footprints']:        
        lyr = arcpy.mapping.ListLayers(mxd, lyrName, df)[0]
        lyr.replaceDataSource(outFolder, "SHAPEFILE_WORKSPACE", "%s%s" % (cntry, lyr.datasetName[3:]))    
        lyr.definitionQuery = ""
        
    focalCountry = arcpy.mapping.ListLayers(mxd, 'Focus Country', df)[0]
    allCountries = arcpy.mapping.ListLayers(mxd, 'Other Countries', df)[0]
    cities = arcpy.mapping.ListLayers(mxd, 'Cities', df)[0]
    focalCountry.definitionQuery = '"ISO3" = \'%s\'' % cntry    
    allCountries.definitionQuery = 'NOT "ISO3" = \'%s\'' % cntry    
    #Define the cities to view on the final map
    getCities = True
    loopCnt = 0
    minPop = 100000
    while getCities:
        cities.definitionQuery = '"ES00POP" > %s AND "ISO3" = \'%s\'' % (minPop, cntry)
        loopCnt += 1
        result = int(arcpy.GetCount_management(cities).getOutput(0))
        tPrint("City population of %s yeilded %s cities" % (minPop, result))
        #If the number of cities is not between 10 and 2, adjust the min population
        if (result < 10 and result > 2) or loopCnt > 9:
            getCities = False
        elif result > 10:
            minPop = minPop * (5./3.)
        else:
            minPop = minPop * (1./3.)
    return( { 'Extent':focalCountry.getExtent(), 'FPPath':arcpy.mapping.ListLayers(mxd, 'Rad Cal Footprint 2010', df)[0].dataSource } )
        
def FPpostProcessing(docsFolder, mapsFolder, gisFolder, outDef, cntry, finalShape, newCities, thresh, popField = "ES00POP", createMaps=False):
    tPrint("All Nighttime Lights Processed, extracting final output")
    #Input maps used to create output maps
    inMxdFootprints = r"S:\GLOBAL\Projects\CityLights\Data\NighttimeLights_Extents2.mxd"
    inMxdCMPD = r"S:\GLOBAL\Projects\CityLights\Data\GlobalNTL_Maps_CMPD2.mxd"
    #Input and output files in final nighttime lights deliverable
    inPPT = r"S:\GLOBAL\Projects\CityLights\Data\Introduction to NTL.pptx"
    inWordDoc = r"S:\GLOBAL\Projects\CityLights\Data\NTL technical documentation.docx"
    outPPT = os.path.join(docsFolder, "%s Introduction to NTL.pptx" % cntry)
    outDoc = os.path.join(docsFolder, "NTL technical documentation.docx")
    if createMaps:
        outMxd = os.path.join(mapsFolder, "%s NighttimeLights_Extents.mxd" % cntry)
        outMxd2 = os.path.join(mapsFolder, "%s NighttimeLights_CMPD.mxd" % cntry)
        outCmpd = os.path.join(mapsFolder, "%s CmpdGrowth_1996-2010_Clipped_Extents.png" % cntry)
        
    #write the output from the master footprints and the cities database to an output xls    
    workbook = openpyxl.Workbook()
    workbook = writeWorkbookDictTitlePage(workbook)
    workbook.save(outDef)

    #rawStyle = workbook.add_format({'font_color':'#006153', 'bg_color':'#C6EFCE'})
    cityIgnoreFields = ["Nothing","nearT0FID","nearT1FID","t0Dist","t1Dist","NEAR_FID","NEAR_DIST","GRIDCODE","URBORRUR",
                        "YEAR","INSGRUSED","CONTINENT","UNREGION","COUNTRY","UNSD","ISO3","SCHNMTYPE","SRCTYP","COORDSRCE",
                        "DATSRC","LOCNDATSRC","NOTES","ADMNM2","TYPE","LATLONGID","SCHADMNM","NUTS2","NUTS_3","POINT_X","POINT_Y"]
    orderFields = [["ExtentName",1],["ExtTypeT0",2],["CtyCntT0",3],["ExtTypeT1",4],["CtyCntT1",5],["Status",6]]
    #Add a unique field to the NTL Extents Shapefile
    arcpy.AddField_management(finalShape, "unqID", "TEXT")
    arcpy.CalculateField_management(finalShape, "unqID", "'%s_%%s' %% !FID!" % cntry, 'PYTHON_9.3')
    
    writeShapefileXLS(outDef, finalShape, "Extents", "pop", "D", ["GRIDCODE"], orderFields)
    writeShapefileXLS(outDef, newCities, "Cities", popField, "D", cityIgnoreFields, orderFields)
    workbook = openpyxl.load_workbook(outDef)
    workbook = writeWorkbookCharts(workbook)
    workbook.save(outDef)

    if not os.path.exists(outPPT):
        shutil.copyfile(inPPT, outPPT)    
    if not os.path.exists(outDoc):
        shutil.copyfile(inWordDoc, outDoc)
    if createMaps:
        #Create the Nighttime Lights Footprint Map    
        mxd = arcpy.mapping.MapDocument(inMxdFootprints)
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        for df in arcpy.mapping.ListDataFrames(mxd):
            outExtent = fixMapConnections(mxd, df, cntry, gisFolder)
            df.extent = outExtent['Extent']           
        #Export final footprint image
        arcpy.mapping.ExportToPNG(mxd, os.path.join(mapsFolder, "%s Nighttime Lights Extents.png" % cntry))
        mxd.saveACopy(outMxd)
        
        #Create maps of compound growth rate
        mxd2 = arcpy.mapping.MapDocument(inMxdCMPD)
        df = arcpy.mapping.ListDataFrames(mxd2)[0]
        allLayers = arcpy.mapping.ListLayers(mxd2, "", df) 
        focalCountry = allLayers[0]
        allCountries = allLayers[1]
        focalCountry.definitionQuery = '"ISO3" = \'%s\'' % cntry
        df.extent = outExtent['Extent']
        allCountries.definitionQuery = 'NOT "ISO3" = \'%s\'' % cntry
        #Create Compound Growth Rate Map for everything
        cmpdLayer = allLayers[6]
        cmpdMask = allLayers[5] #This layer masks out the cmpd growth map to just the cities
        cmpdLayer.visible = True
        cmpdMask.visible = False
        arcpy.mapping.ExportToPNG(mxd2, os.path.join(mapsFolder, "%s %s" % (cntry, "CompoundGrowthRate_ALL.png")))
        #Create CMPD map for just the city footprints
        cmpdMask.visible = True
        cmpdMask.symbology.classBreakValues = [0, thresh, 7000]
        arcpy.RefreshActiveView()
        arcpy.mapping.ExportToPNG(mxd2, os.path.join(mapsFolder, "%s %s" % (cntry, "CompoundGrowthRate_Extents.png")))

    #Attribute Metadata for fpT0, fpT1, fpResults["DISFile"], fpResults["APPFile"], finalShape, newCities
    GOSTRocks.metadataManagement.setMetadata("%s/%s_%s.xml" % (gisFolder, cntry, "radCalT0.shp"), {'Title': "Nighttime Lights Extents 1996", 
        'Summary':"Urban Extents extracted from Nighttime Lights", 
        'Description': "Urban extents in %s were extracted from Nighttime Lights data (http://ngdc.noaa.gov/eog/dmsp/download_radcal.html) in %s  \
by applying a threshold of %s, calculated through a comparison with ESA globcover urban maps." % (cntry, 1996, thresh),
        'Keywords': "GLOBAL, Urban, NTL, Nighttime Lights, %s" % cntry, 
        'Limitations': "For free, public use",
        'Credits': 'WBG, Benjamin P. Stewart'})   
    GOSTRocks.metadataManagement.setMetadata("%s/%s_%s.xml" % (gisFolder, cntry, "radCalT1.shp"), {'Title': "Nighttime Lights Extents 2010", 
        'Summary':"Urban Extents extracted from Nighttime Lights", 
        'Description': "Urban extents in %s were extracted from Nighttime Lights data (http://ngdc.noaa.gov/eog/dmsp/download_radcal.html) in %s  \
by applying a threshold of %s, calculated through a comparison with ESA globcover urban maps." % (cntry, 2010, thresh), 
        'Keywords': "GLOBAL, Urban, NTL, Nighttime Lights, %s" % cntry, 
        'Limitations': "For free, public use",
        'Credits': 'WBG, Benjamin P. Stewart'})       
    GOSTRocks.metadataManagement.setMetadata("%s/%s_%s.xml" % (gisFolder, cntry, "masterExtents.shp"), {'Title': "Nighttime Lights Extents", 
        'Summary':"Urban Extents extracted from Nighttime Lights", 
        'Description': "Urban extents in %s were extracted from Nighttime Lights data (http://ngdc.noaa.gov/eog/dmsp/download_radcal.html) in %s  \
by applying a threshold of %s, calculated through a comparison with ESA globcover urban maps. The same process was performed on the 2010 data. The final\
extents were compared to an inventory of cities to determine footprint names, populations, and status. This file contains all the extents that were \
identified in both time periods, as well as the disappearing and appearing extents. For all features, a number of attributes were calculated taht \
describe the changing nature of the city's brightness." % (cntry, 1996, thresh), 
        'Keywords': "GLOBAL, Urban, NTL, Nighttime Lights, %s" % cntry, 
        'Limitations': "For free, public use",
        'Credits': 'WBG, Benjamin P. Stewart'})

def writeWorkbookDictTitlePage(workbook):
    #Write a Title Page
    #titleCenter = workbook.add_format({'font_color':'white', 'bg_color':'#4D1434', 'bold':True, 'font_size':16, 'center_across':True})    
    #titleFormat = workbook.add_format({'font_color':'white', 'bg_color':'#4D1434'})           
    
    titlePage = workbook.active
    titlePage.title = "Title Page"
    #titlePage.set_column(0,0,130)
    titlePage.cell(row=1,column=1).value = "Global Nighttime Lights Urban Extents and Growth Product - Alpha Release"
    titlePage.cell(row=2,column=1).value = ""   
    titlePage.cell(row=3,column=1).value = "South Asia Urban Team, GSURR: Mark Roberts and Mihir Prakash"
    titlePage.cell(row=4,column=1).value = "In technical collaboration with ITSOP: Benjmain P Stewart and Katie McWilliams"
    titlePage.cell(row=5,column=1).value = ""
    titlePage.cell(row=6,column=1).value = "This Spreadsheet contains 5 pages"
    titlePage.cell(row=7,column=1).value = "1. Title Page"
    titlePage.cell(row=8,column=1).value = "2. Data Dictionary - Describes the column definitions in sheets 3 and 4"
    titlePage.cell(row=9,column=1).value = "3. Urban Extents - Summary of urban areas extracted from radiance calibrated nighttime lights"
    titlePage.cell(row=10,column=1).value = "4. Cities Inventory - Cities inventory with descriptive nighttime lights attributes (default Cities are from GRUMP 2000 populated places inventory)"
    titlePage.cell(row=11,column=1).value = "5. Charts - Sample charts showing change in brightness in 10 largest urban extents."
    
    #Write a Data Dictionary
    dataDict = workbook.create_sheet("Data Dictionary")
    dataDict.cell(column=2,row=1).value = "Extents and Cities Column Definitions"
    dataDict.cell(column=1,row=1).value = ""
    #dataDict.set_column(1,1,1200)
    rIdx = 1
    for rVal in ["FID","ExtentName","ExtTypeT0","CtyCntT0","ExtTypeT1","CtyCntT1","Status","gAreaKM","pop","RC1996_T0","RC2010_T1","ntlChange","ntlChgCorr","intensive","extensive","extenCorr","areaChg","raw1992","rc1996"]:
        rIdx = rIdx + 1
        dataDict.cell(row=rIdx,column=1).value = rVal
    
    rIdx = 1
    for rVal in ["GIS Code for joining to Master Extents shapefiles",
        "Name of extent if found",
        "Type of extent in 1996 (Agglomeration, Stand-alone extent, -1 found but without an intersecting city, blank indicates there was no comparable extent)",
        "Number of cities in urban extent in 1996 (only greater than 1 for agglomerations)",
        "Type of extent in 2010 (Agglomeration, Stand-alone urban, -1 found but without an intersecting city, blank indicates there was no comparable extent)",
        "Number of cities in urban extent in 2010 (only greater than 1 for agglomerations)",
        "Status of urban extent - Missed (does not intersect any cities) Found (intersects cities in both time periods) Disappear (intersects city only in 1996) Appear (intersect city only in 2010)",
        "Area of urban extent in 2010 in km2",
        "Population tabulated from GRUMP data  (~2000)",
        "Total brightness of 1996 radiance calibrated NTL for urban extent in T0",
        "Total brightness of 2010 radiance calibrated NTL for urban extent in T1",
        "Change of brightness in urban extent (RC2010_T1 - RC1996_T0)",
        "Change in brightness of urban extent corrected for existing light in expanded area (RC2010_T1 - RC1996_T0 - (RC1996_T1 - RC1996_T0))",
        "Change of brightness in urban core (RC2010_T0 - RC1996_T0)",
        "Change of brightness in newly created area of urban extent (RC2010_T1 - RC2010_T0)",
        "Change of brightness in newly created area of urban extent corrected for existing light in expanded area (RC2010_T1 - RC2010_T0 - (RC1996_T1 - RC1996_T0))",
        "Change in area from 1996 to 2010 in km2",
        "Total brightness in raw NTL in 1992 for T1 extent",
        "Total brightness in radiance calibrated NTL in 1992 for T1 extent"]:
        rIdx = rIdx + 1
        dataDict.cell(row=rIdx,column=2).value = rVal 
    return workbook
    
def writeWorkbookCharts(workbook):
    chartPage = workbook.create_sheet("Charts")
    extentChart = openpyxl.chart.BarChart()
    extentChart.type = "col"
    extentChart.style = 10
    extentChart.title = 'Change in Extent Brightness'
    extentChart.x_axis.title = "10 largest cities"
    extentChart.y_axis.title = 'Total change in brightness'
    
    #Add data on NTL Change - Total, Intensive, and Extensive
    dataSheet = workbook['Extents']
    data = openpyxl.chart.Reference(dataSheet, min_col=13, max_col=15, min_row=1, max_row=10)
    cats = openpyxl.chart.Reference(dataSheet, min_col=2, min_row=1, max_row=10)
    extentChart.add_data(data)
    extentChart.set_categories(cats)    
    chartPage.add_chart(extentChart, "A1")
    return workbook
'''
    chartPage.insert_chart('A1', extentChart)
    
    #Add line chart for each of the cities
    #Add data on NTL Change - Total, Intensive, and Extensive
    lineChart = workbook.add_chart({'type':'line'})
    for idx in range(2,11):
        lineChart.add_series({'name':'=Extents!$B%s' % idx,'values':'=Extents!$R$%s:$AL$%s' % (idx,idx),'categories':'=Extents!$R$1:$AL$1'})
    #Set titles
    lineChart.set_title({'name':'Total Extent Brightness (Raw NTL)'})
    lineChart.set_x_axis({'name':'Years of Raw NTL Brightness'})
    lineChart.set_y_axis({'name':'Total Brightness','min':0})
    lineChart.set_size({'width':600,'height':600})
    chartPage.insert_chart('K1', lineChart)
    
    #Add line chart for each of the cities
    #Add data on NTL Change - Total, Intensive, and Extensive
    rcChart = workbook.add_chart({'type':'line'})    
    for idx in range(2,11):        
        rcChart.add_series({'name':'=Extents!$B%s' % idx,'values':'=Extents!$AM$%s:$AS$%s' % (idx,idx),'categories':'=Extents!$AM$1:$AS$1'})
    #Set titles
    rcChart.set_title({'name':'Total Extent Brightness (RadCal NTL)'})
    rcChart.set_x_axis({'name':'Year'})
    rcChart.set_y_axis({'name':'Total Brightness','min':0})
    rcChart.set_size({'width':600,'height':600})
    chartPage.insert_chart('U1', rcChart)
'''