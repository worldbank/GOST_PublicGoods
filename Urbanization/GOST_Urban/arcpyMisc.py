###############################################################################
# Geospatial Operational Support Team @ The World Bank Group
# Benjamin P. Stewart
# Purpose: miscellaneous modules 
###############################################################################

import time, sys, os, csv, re, xlsxwriter, math, glob, arcpy, json, inspect
import pandas as pd
from dbfpy import dbf
from GOSTRocks.misc import tPrint

from arcpy.sa import *

arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

def getUrbanParams():
    thisFolder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    inputParameters = "%s/Urban/urbanParameters.json" % thisFolder
    with open(inputParameters, 'rb') as data_file:                 
        jsonParams = json.load(data_file)
    return jsonParams

'''
Returns an arcpy Field object from the selected layer
---inLyr: string path to file
---name: name of desired field
'''
def getFieldObject(inLyr, name):
    fields = arcpy.ListFields(inLyr)
    for f in fields:
        if f.name == name:
            return(f)
'''
Returns the sum of the defined field
'''
def getFieldSum(inShp, fieldName):
    cur = arcpy.da.SearchCursor(inShp, [fieldName])
    total = 0
    for c in cur:
        total = total + c[0]
    return (total)
            
'''
Rename a field in a feature class
'''
def renameField(inShp, oldFieldName, newField):
    #Find old field
    oldField = getFieldObject(inShp, oldFieldName)
    if not oldField:
        return -1    
    #Add new field to shapefile with type matching old field
    tryAddField(inShp, newField, oldField.type)    
    #Calculate old values into new shapefile
    arcpy.CalculateField_management(inShp, newField,"!%s!" % oldFieldName,"PYTHON_9.3")
    #Delete old Field
    arcpy.DeleteField_management(inShp, oldFieldName)
    return(1)

'''
Add a field to a feature class if it doesn't exist
'''    
def tryAddField(inShp, fieldName, fieldType="DOUBLE"):    
    for f in arcpy.ListFields(inShp):
        if f.name == fieldName:
            return(-1)
    arcpy.AddField_management(inShp, fieldName, fieldType)

'''
return a list of unique values in an input table
'''    
def uniqueValues(table, field):
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})

'''
Count of the number of features within a table
'''
def getFeatureCount(inFC, lyrName="EleshNorn"):
    arcpy.MakeTableView_management(inFC, lyrName)
    return(int(arcpy.GetCount_management(lyrName).getOutput(0)))

    
'''Open all the dbf files in dbf list and extract all relevant values
---dbfList: list of dbf files to process, these should be from zonal statistics in most cases
---idField: Field that identifies the individual zonal features
---sumVals: List of equal length to dbf list describing the fields to extract from dbf. Use -1 to get just the sums
-RETURNS: A data dictionary of idFields with
'''
def summarizeDbf(dbfList, idField, sumVals=-1):
    results = {}
    expectedLength = 0
    if sumVals == -1:
        sumVals = [["SUM"]] * len(dbfList)
    else:
        sumVals = sumVals #[sumVals] * len(dbfList)
    for dFileIdx in range(0, len(dbfList)):
        dFile = dbfList[dFileIdx]
        dbfData = dbf.Dbf(dFile)
        for idx in range(0, len(dbfData)):
            curResults = []
            for cVal in sumVals[dFileIdx]:
                curResults.append(dbfData[idx][cVal])
            if dFile == dbfList[0]:
                results[dbfData[idx][idField]] = curResults
            else:                
                try: #For the 2nd dbf and onward, append the results to the existing data
                    #Check to make sure the current data are of the proper length, if not, pad with zeroes
                    curRes = results[dbfData[idx][idField]]
                    if len(curRes) != expectedLength:          
                        results[dbfData[idx][idField]].extend([0] * (expectedLength - len(curRes)))
                    results[dbfData[idx][idField]].extend(curResults)                    
                except:
                    #If the key is not in the output, add a bunch of zeroes at the start, then append
                    results[dbfData[idx][idField]] = [0] * expectedLength
                    results[dbfData[idx][idField]].extend(curResults)
        expectedLength = expectedLength + len(sumVals[dFileIdx])
    return(results)   

'''
Write a dictionary to an output file (meant to be used with function above
---inD: input dictionary
---outFile: output csv file
---headers: column titles, written as first row of the output file
'''
def writeDict(inD, outFile, headers):            
    outputFile = open(outFile, 'wb')    
    outputFile.write(','.join(str(h) for h in headers) + "\n")
    for key, value in inD.items():
        outputFile.write("%s,%s\n" %(key, ','.join(str(f) for f in value)))
    outputFile.close()

'''
Extract the file names from a VRT stack
'''
def getVRTTitles(inVRT):
    inD = open(inVRT, 'r')
    allTitles = []
    for line in inD:
        if 'SourceFilename' in line:
            inType = re.search('VRT="1">(.*)/', line).group(1).split("/")[0]
            inScale = re.search('_sc(.*)_fea', line).group(1)
            inFeat = re.search('_fea(.*).tif', line).group(1)           
            inT = "%s_%s_%s" % (inType, inScale, inFeat)
            allTitles.append("%s_%s" % (inT, 'SUM'))
            allTitles.append("%s_%s" % (inT, 'MEAN'))
            allTitles.append("%s_%s" % (inT, 'STD'))
    return(allTitles)
       
'''
Get the unique values from an input raster dataset
'''    
def getUniqueValues(fR):
    try:
        fR = Int(Raster(fR))
        arcpy.BuildRasterAttributeTable_management(fR, "NONE")
        unqVals = []
        for r in arcpy.SearchCursor(fR):
            unqVals.append(r.getValue('Value'))
        return(unqVals)
    except:
        raise(TypeError("The input categorical raster cannot have a raster table", fR))

'''
get a list of nightime light raster files
'''
def getNTLFiles(rawFolder = r"S:\GLOBAL\NightLights",radCalFolder = r"S:\GLOBAL\NightLights\rad_cal"):      
    ntlFiles = []
    for yr in range(1992, 2014):
        curFolder = os.path.join(rawFolder, str(yr))
        inFiles = os.listdir(curFolder)
        for f in inFiles: 
            if re.search("ElvidgeCorrected_gt3.tif$", f):
                ntlFiles.append(os.path.join(curFolder, f))
    radCalFiles = os.listdir(radCalFolder)    
    for f in radCalFiles:
        if re.search("_vis_Corrected.tif$", f):
            ntlFiles.append(os.path.join(radCalFolder,f))
    return(ntlFiles)

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
            
'''
Run zonal stats on GHSL based on the input inShp
'''    
def summarizeGHSL(inShp, idField, threshold=50, part=3):
    inFiles = [r"S:\GLOBAL\GlobalUrbanSettlement\GHSL_Landsat_1975_aggregated_300m\Part%s.tif" % part,
                r"S:\GLOBAL\GlobalUrbanSettlement\GHSL_Landsat_1990_aggregated_300m\Part%s.tif" % part,
                r"S:\GLOBAL\GlobalUrbanSettlement\GHSL_Landsat_2000_aggregated_300m\Part%s.tif" % part,
                r"S:\GLOBAL\GlobalUrbanSettlement\GHSL_Landsat_2014_aggregated_300m\Part%s.tif" % part]
    summarizeVals = []
    filesToDelete = []    
    for fIdx in range(0, len(inFiles)):
        tempRast = "C:/Temp/aaaTemp%s.tif" % fIdx
        outTable = "C:/Temp/ghslTable%s.dbf" % fIdx
        inFile = inFiles[fIdx]
        #Clip the input image to the input shapefile
        arcpy.Clip_management(inFile, "#", tempRast, inShp, '', "ClippingGeometry")
        #Convert the input data to a binary raster based on the input threshold
        tRast = Raster(tempRast)
        tRast = tRast > threshold
        #Run zonal statistics on binary data
        ZonalStatisticsAsTable(inShp, idField, tRast, outTable)
        summarizeVals.append(["SUM"])
        filesToDelete.append(outTable)
    if idField in ['FID', 'OID']:
        idField = "%s_" % idField
    res = summarizeDbf(filesToDelete, idField, summarizeVals)
    return({"Results":res, "Titles":["ID","GHSL1975","GHSL1990","GHSL2000","GHSL2014"]})    

'''
Run zonal statistics on a series of base global datasets and return results
    the inFiles can be replaced by your own versions if desired
'''
def summarizeGlobalData(inShp, idField, otherFiles=-1, type=['N','N','N','N','N','N','N','C'], clip=True, tempFolder="C:/Temp"):
    inFiles = [r"S:\GLOBAL\NightLights\1992\F101992.v4b_web.stable_lights.avg_vis_ElvidgeCorrected.tif",
               r"S:\GLOBAL\NightLights\2012\F182012.v4c_web.stable_lights.avg_vis_ElvidgeCorrected.tif",
               r"S:\GLOBAL\NightLights\rad_cal\F1996_0316-19970212_rad_v4.avg_vis_Corrected.tif",
               r"S:\GLOBAL\NightLights\rad_cal\F2010_0111-20101209_rad_v4.avg_vis_Corrected.tif",
               r"S:\GLOBAL\GDP\UNEP\GDP.tif",
               r"S:\GLOBAL\Elevation\1km\elevation1",
               r"S:\GLOBAL\POPULATION\Landscan2012\ArcGIS\Population\lspop2012",
               r"S:\GLOBAL\Landcover\Globcover\GLOBCOVER_L4_200901_200912_V2.3.tif"]
    inTitles = [idField, "NTL2012", "RC1996", "RC2010", "GDP_UNEP", "Elev", "L11","L14","L20","L30","L40",
                "L50","L60","L70","L90","L100","L110","L120","L130","L140","L150","L160","L170","L180",
                "L190","L200","L210","L220","L230","NTL1992"]
    if not otherFiles == -1 and not otherFiles == "NTL":
        inFiles = otherFiles
        inTitles = []
    if otherFiles == "NTL":
        inFiles = getNTLFiles()
        type = ["N"] * len(inFiles)
    inTitles = [idField]
          
    #Define the values to extract from each DBF
    summarizeVals = []
    filesToDelete = []
    tempVals = []
    fCount = 0
    for inFidx in range(0, len(inFiles)):
        tempRast = os.path.join(tempFolder, "aaaTemp%s.tif" % inFidx)
        tempVals.append(tempRast)
        inFile = inFiles[inFidx]        
        fName = os.path.basename(inFile)        
        curType = type[inFidx]
        if not clip:
            tempRast = inFile
        else:
            arcpy.Clip_management(inFile, "#", tempRast, inShp, '', "ClippingGeometry")
        inR = Raster(tempRast)
        if curType == 'C':
            #Processing categorical data is slightly different             
            unqVals = getUniqueValues(inFile)
            print(unqVals)
            for unq in unqVals:
                summarizeVals.append(["SUM"])  #All the values from categorical datasets need SUM
                #For the current unique value, create a binary raster
                curR = inR == unq
                outTable = os.path.join(tempFolder, "%s_%s.dbf" % (fName.replace(".", "").replace(".tif", ""), unq))
                inTitles.append("%s_%s" %(fName, unq))
                ZonalStatisticsAsTable(inShp, idField, curR, outTable, 'DATA')
                filesToDelete.append(outTable)     
        else:
            summarizeVals.append(["SUM","MEAN","STD"])
            inTitles.append("%s_SUM" % fName)
            inTitles.append("%s_MEAN" % fName)
            inTitles.append("%s_STD" % fName)
            
            outTable = os.path.join(tempFolder, "%s_%s.dbf" % (fName.replace(".", "").replace(".tif", "").replace("-", "_"), fCount))
            fCount += 1
            if not os.path.exists(outTable):
                ZonalStatisticsAsTable(inShp, idField, Raster(tempRast), outTable)
            filesToDelete.append(outTable)
            tPrint("Finshed Processing " + os.path.basename(inFile).replace(".tif", ""))
        
    if idField in ['FID', 'OID']:
        idField = "%s_" % idField
    res = summarizeDbf(filesToDelete, idField, summarizeVals)

    for f in filesToDelete:
        arcpy.Delete_management(f)    
    for f in tempVals:
        arcpy.Delete_management(f)
    return({"Results":res, "Titles":inTitles})
    
'''Convert raster values above the defined threshold into polygons
---tempRaster: input raster to be converted
---threshold: cutoff for raster polygons. Above is polygon, below is not
---outPolygons - output vector data file
'''
def extractFootprints(tempRaster, threshold, outPolygons, selectTerm=' "GRIDCODE" > 0 '):
    tempShp = "C:/Temp/cntry.shp"
    
    tRaster = Raster(tempRaster)
    tRaster = tRaster > int(float(threshold))

    #Convert to vector and select out the non zero values
    arcpy.RasterToPolygon_conversion(tRaster, tempShp, 'NO_SIMPLIFY', '#')
    arcpy.Select_analysis(tempShp, outPolygons, selectTerm)
    #Add area field and calculate area in square kilometres
    arcpy.AddField_management(outPolygons, "gAreaKM", "FLOAT")
    arcpy.CalculateField_management(outPolygons, "gAreaKM", "!SHAPE.AREA@SQUAREKILOMETERS!", "PYTHON")
    arcpy.Delete_management(tempShp)

def createMapFromColumns(mxd, layerName, outputFolder, columns="all"):
    '''
    Create a series of maps, one for each column in the layer defined by layerName
    mxd - Input mapDocument object
    layername - name of layer to loop through
    outputFolder - place to create maps

    columns - array of columns to create map for
    
    ###LOOK here for a QGIS example
    http://www.qgis.nl/2013/08/13/python-script-to-generate-series-of-maps/?lang=en
    '''
    df = arcpy.mapping.ListDataFrames(mxd)[0]
    #Get refernece to focal layer
    layers = arcpy.mapping.ListLayers(mxd, '', df)
    for l in layers:
        print l.name
        if l.name == layerName:
            featureLayer = l
    #Get list of columns in focal layer
    print featureLayer
    fieldNames = [f.name for f in arcpy.ListFields(featureLayer)]
    if columns != 'all':
        fieldNames = columns
    
    for f in fieldNames:
        featureLayer.sympology.valueField = f
        featureLayer.sympology.numClasses = 8
        featureLayer.symbology.reclassify()
        mxd.save()
        time.sleep(1)
        arcpy.mapping.ExportToPNG(mxd, os.path.join(outputFolder, "%s_%s.png" % (layerName, f)))
        
        
    
    
#Create a map from a mxd
###mxd - Path to input MXD
###outputImage - the output Image to create
###visibility - binary variables to define the layers as on
#   or off. If the visibility list is shorter than the number of layers, the remaining
#   follow the definition from the original Map
###statusDefinition - Dictionary containing status defintions, key in dictionary is layername
###extent - set the extent of the output map as an extent feature
###scale - set the scale of the map as a number
def createMapFromMXD(mxd, outputImage, visibility = [], statusDefinition={}, extent="", scale=""):    
    df = arcpy.mapping.ListDataFrames(mxd)[0]
    if extent != "":
        df.extent = extent
    if scale != "":
        df.scale = scale
   
    if len(visibility) > 0:
        layers = arcpy.mapping.ListLayers(mxd, '', df)
        for lIdx in range(0, len(visibility)):
            layers[lIdx].visible = visibility[lIdx] 
            # print ("%s - %s" % (layers[lIdx].name, layers[lIdx].visible))
    arcpy.mapping.ExportToPNG(mxd, outputImage)

#Create a map from a mxd for each feature in a defined class
###mxd - Path to input MXD
###outputImage - the output Image to create
###visibility - binary variables to define the layers as on
#   or off. If the visibility list is shorter than the number of layers, the remaining
#   follow the definition from the original Map
###statusDefinition - Dictionary containing status defintions, key in dictionary is layername
###extent - set the extent of the output map as an extent feature
###scale - set the scale of the map as a number
def createMapFromMXD_loop(mxd, outputImage, zoomLayer, layerName, visibility = [], statusDefinition={}, extent="", scale=""):
    df = arcpy.mapping.ListDataFrames(mxd)[0]
    if extent != "":
        df.extent = extent
    if scale != "":
        df.scale = scale
   
    if len(visibility) > 0:
        for lIdx in range(0, len(visibility)):
            layers[lIdx].visible = visibility[lIdx]
    #Loop through the features in one of the feature classes
    layers = arcpy.mapping.ListLayers(mxd)
    for l in layers:        
        #If there is a staus definition, apply it here
        if l.name in statusDefinition.keys():
            l.definitionQuery = statusDefinition[l.name]
            print ("Definition Query set for %s" % l.name)
        if l.name == zoomLayer:
            loopLyr = l
            l.visible = True
            
    with arcpy.da.SearchCursor(loopLyr, ["OID@", layerName]) as cur:
        for feat in cur:
            #try:
            print feat[1]
            arcpy.SelectLayerByAttribute_management(loopLyr, "NEW_SELECTION", '"FID" = %s' % feat[0])
            #df.zoomToSelectedFeatures()
            df.extent = loopLyr.getSelectedExtent()
            #time.sleep(5)
            arcpy.SelectLayerByAttribute_management(loopLyr, "CLEAR_SELECTION", '"FID" = %s' % feat[0])
            outputMap = outputImage.replace(".png", "_%s.png" % feat[1])                  
            arcpy.mapping.ExportToPNG(mxd, outputMap)
            #except:
            #    print("Something failed with feature %s" % feat[1])

def getLargestPointsInPolygons(inPolys, inPoints, weightField, polyIdField, ptsIdField):
    '''
    For each Polygon in inPolys, return the point with the largest value in weightField
    inPolys - polygons within which to search for points
    inPoints - points to search for within the inPolys
    weightField - gravity field to search for in the inPoints
    idField - id field for inPoints
    RETURNS - dataframe of cities with the following columns ptIDField, adminIDField, weightField, lat, long
    '''
    outputDF = pd.DataFrame(columns=["ptID", "polyID", "pop", "lat", "lng"])
    #outputDF = []
    pointLayer = "pyLyr"
    arcpy.MakeFeatureLayer_management(inPoints, pointLayer)
    featCnt = 0
    with arcpy.da.SearchCursor(inPolys, ["SHAPE@XY", "OID@", polyIdField]) as adminCur:
        for feature in adminCur:            
            featCnt = featCnt + 1
            #Create a feature layer for the first feature from the OID and select intersecting points
            adminLyr = "FUBAR_%s" % feature[2]
            arcpy.MakeFeatureLayer_management(inPolys, adminLyr, '"FID" = %s' % feature[1])
            arcpy.SelectLayerByLocation_management(pointLayer, "INTERSECT", adminLyr)
            count = int(arcpy.GetCount_management(pointLayer)[0])
            if featCnt % 10 == 0:
                tPrint("Polygon %s intersects %s points" %(feature[2], count))
            maxWeight = -1
            #Search through all point features that intersect the current admin shape
            with arcpy.da.SearchCursor(pointLayer, ["SHAPE@XY", "OID@", weightField, ptsIdField]) as ptCursor:
                for ptFeat in ptCursor:
                    if ptFeat[2] > maxWeight:
                        maxWeight = ptFeat[2]
                        outputFeature = ptFeat
            #Append the selected feature to the output
            if maxWeight > 0:                
                curDF = pd.DataFrame([[outputFeature[3], feature[2], outputFeature[2], outputFeature[0][0], outputFeature[0][1]]], 
                                     columns=["ptID", "polyID", "pop", "lat", "lng"])         
            else: #If no points intersect area, use the admin's centroid
                curDF = pd.DataFrame([["CENTROID", feature[2], 0, feature[0][0], feature[0][1]]], 
                                     columns=["ptID", "polyID", "pop", "lat", "lng"]) 
            outputDF = outputDF.append(curDF)
            if featCnt > 20:
                return outputDF
                
    #return pd.concat(outputDF, axis=2)
    return outputDF
    
def pointsInPolygons(inPolys, inShapes, fieldName, attachFields, inFields=["SCHNM","MAX","ES00POP"], 
    outField="CityName", clipPts=False, clipFile="D:/Temp/pts.shp", whereClause=""):
    '''    
    Calculate the intersecting inShapes that intersect inPolys, optional fields allow for intersecting values to be calculated
    inPolys: the input shapes to be updated with intersection numbers
    inShapes: comparison shapefile (works with any geometry)
    attachFields: Boolean determining if fields should be added to inPolys

    ##BEN## - Default is to attach cityname from GRUMP city shapefile by maximum 2000 population estimate
    inFields: list of input fields in inShapes that should be joined to inPolys along with joining procedure
    -- MAX/MIN indicates the value of the highest value in the matching field is attached
    -- SUM indicates the intersecting field values should be summed
    outFields: output fields to be added to match inFields
    clipPts: Performance variable - should be set to true if the inShapes is lareg geographically (ie-global)
             and the inPolys are small (national, regional, etc)
    '''      
    tempLyr = "FUBAR_Lyr"
    citiesLyr = "CITIES"
    cursorFields = ["OID@", "SHAPE@", fieldName]
    tryAddField(inPolys, fieldName) 
        
    #[OPTIONAL] ClipPoints
    if clipPts:
        arcpy.Clip_analysis(inShapes, inPolys, clipFile)
        inShapes = clipFile
    arcpy.MakeFeatureLayer_management(inShapes, citiesLyr)
    
    #[OPTIONAL] Add output fields
    if attachFields:        
        inField = getFieldObject(inShapes, inFields[0])
        tryAddField(inPolys, outField, inField.type)
        cursorFields = ["OID@", "SHAPE@", fieldName, outField]

    #Run update cursor on inPolys
    with arcpy.da.UpdateCursor(inPolys, cursorFields,"") as cur:
        for row in cur:        
            #Make feature layer of the current feature
            arcpy.MakeFeatureLayer_management(inPolys, tempLyr, '"FID" = %s' % row[0])
            arcpy.SelectLayerByLocation_management(citiesLyr, 'INTERSECT', tempLyr)
            count = int(arcpy.GetCount_management(citiesLyr)[0])
            row[2] = count
            #If fields from the intersection are to be attached, do it here
            if attachFields:
                minVal = 1000000000
                maxVal = 0
                sumVal = 0 
                with arcpy.da.UpdateCursor(citiesLyr, [inFields[0], inFields[2]]) as ptCur:
                    for ptRow in ptCur:                                          
                        if ptRow[1] > maxVal and inFields[1] == "MAX":
                            maxVal = ptRow[1] 
                            row[3] = ptRow[0]
                        if ptRow[1] < minVal and inFields[1] == "MIN":
                            minVal = ptRow[1] 
                            row[3] = ptRow[0]
                        if inFields[1] == "SUM":
                            sumVal = sumVal + ptRow[1] 
                            row[3] = sumVal
            cur.updateRow(row)
            print row
    
    

def kappa(inRaster1, inRaster2, tempFolder=r"D:\Temp"):    
    #Open raster, convert to binary urban and non-urban
    inR1 = Raster(inRaster1) > 0
    inR1 = Con(IsNull(inR1), 0, 1)
    #Open raster, convert to binary urban and non-urban
    inR2 = (Raster(inRaster2) > 0) * 10
    inR2 = Con(IsNull(inR2), 0, 10)
    
    inK = inR1 + inR2            
    inK.save(os.path.join(tempFolder, "tempR1.tif"))

    arcpy.BuildRasterAttributeTable_management(inK, "NONE")
    
    unqVals = []
    for r in arcpy.SearchCursor(inK):
        unqVals.append([r.getValue('Value'), r.getValue('Count')])
    calculateKappa(unqVals)
    
#Will calculate a single kapp statistic based on a sumarized DBF as created from the above analysis
def calculateKappa(kVals):
    for row in kVals:
        if row[0] == '0':
            aRural = float(row[1])
        if row[0] == '1':
            r1Urban = float(row[1])
        if row[0] == '10':
            r2Urban = float(row[1])
        if row[0] == '11':
            aUrban = float(row[1])
    totalCells = aRural + r1Urban + r2Urban + aUrban
    r1Urban = (r1Urban + aUrban) / totalCells
    r1Rural = (aRural + r2Urban) / totalCells
    r2Urban = (r2Urban + aUrban) / totalCells
    r2Rural = (aRural + r1Urban) / totalCells
    
    pAurban = aUrban / totalCells
    pcAurban = r1Urban * r2Urban
    
    pArural = aRural / totalCells
    pcArural = r1Rural * r2Rural
    
    pAgreeAll = pAurban + pArural
    pcAgreeAll = pcAurban + pcArural
    
    kappa = (pAgreeAll - pcAgreeAll) / (1-pcAgreeAll)
    return(kappa)