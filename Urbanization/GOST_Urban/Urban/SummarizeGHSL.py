################################################################################
# Summarize GHSL
# Benjamin P Stewart, Aug 2015
# Purpose: Run zonal statistics (amongst other things) on GHSL data
################################################################################
import sys, arcpy, os, glob, shutil, xlsxwriter, inspect
import numpy as np
cmd_folder = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
from GOSTRocks.misc import *
from GOSTRocks.arcpyMisc import *
from GOSTRocks.xlsxStuff import *

from arcpy.sa import *

def createYearlyLayers(ghslRaster, outputFolder):
    masterR = Raster(ghslRaster)
    returnList = []
    for yr in [6,5,4,3]:
        curR = masterR >= yr
        outRaster = "%s/ghsl_%s.tif" % (outputFolder, yr)
        returnList.append(outRaster)
        if not arcpy.Exists(outRaster):
            curR.save(outRaster)
    return(returnList)
        

def createSingleLayer(ghslTiles, outFolder, outName):
    arcpy.MosaicToNewRaster_management(";".join(ghslTiles), outFolder, outName, number_of_bands=1)

def summarizeGUF_toShp(inShape, gufTiles, gufInputFolder, tempFolder):
    lyrName = 'gufLyr'
    newFields = ['totalBuilt']
    selFields = ['FID'] + newFields
    fcCount = getFeatureCount(inShape)
    for f in newFields:
        tryAddField(inShape, f, "DOUBLE")
    with arcpy.da.UpdateCursor(inShape, selFields) as cursor:
        for row in cursor:
            featIdx = row[0]
            if featIdx % 10 == 0:
                tPrint("Feature %s of %s" % (featIdx, fcCount))
            arcpy.MakeFeatureLayer_management(inShape, lyrName, '"FID" = %s' % featIdx)
            #Get intersecting GHSL files
            curTiles = selectGHSLtiles(lyrName, gufTiles, 'FileName', False)
            totalSum = 0
            for curFile in curTiles: 
                #for each of the intersecting tiles, clip the tils amd sum the results
                curGHS = os.path.join(gufInputFolder, curFile)
                tempClipped = os.path.join(tempFolder, curFile.replace(".tif", "%s.tif" % featIdx))
                arcpy.Clip_management(curGHS, "#", tempClipped, lyrName, "0", "ClippingGeometry")
                gufArr = arcpy.RasterToNumPyArray(tempClipped, nodata_to_value = 0)
                totalSum += (gufArr > 0).sum()
            row[1] = totalSum
            cursor.updateRow(row)
            
def summarizeGHSL_toShp(inShape, ghslFolder, datamaskFolder, ghslOutputFolder, tempFolder, ghslOutline):
    lyrName = "featLyr"
    newFields = ["totalNB","totalWater","built2014","built2000", "built1990", "built1975"]
    selFields = ["FID"] + newFields
    fcCount = getFeatureCount(inShape)
    for f in newFields:
        tryAddField(inShape, f, "DOUBLE")
    #create update cursor
    allOutFiles = []
    with arcpy.da.UpdateCursor(inShape, selFields) as cursor:
        for row in cursor:
            featIdx = row[0]
            if featIdx % 10 == 0:
                tPrint("Feature %s of %s" % (featIdx, fcCount))
            outFile = r"%s/SumGHSL_%s.csv" % (ghslOutputFolder, featIdx)
            allOutFiles.append(outFile)
            arcpy.MakeFeatureLayer_management(inShape, lyrName, '"FID" = %s' % featIdx)
            #Get intersecting GHSL files
            ghslFiles = selectGHSLtiles(lyrName, ghslOutline, 'Path', False)
            
            for curFile in ghslFiles:               
                builtFile = os.path.join(ghslFolder, curFile)
                openGHSL(builtFile)
                maskFile = os.path.join(datamaskFolder, curFile.replace("MT","DATAMASK").replace('mt', 'datamask'))
                #Create Temp Clipped Files
                tBuilt = os.path.join(tempFolder, os.path.basename(curFile))
                tDataMask = os.path.join(tempFolder, os.path.basename(curFile.replace("mt","datamask")))
                arcpy.Clip_management(builtFile, "#", tBuilt, lyrName, "0", "ClippingGeometry")
                arcpy.Clip_management(maskFile, "#", tDataMask, lyrName, "0", "ClippingGeometry")

                curDict = summarizeMaskByBuilt(tBuilt, tDataMask)
                if curFile == ghslFiles[0]:
                    finalDict = curDict
                else:
                    for kIdx in finalDict:
                        finalDict[kIdx] = [x+y for x,y in zip(finalDict[kIdx], curDict[kIdx])]
            writeDict(finalDict, outFile, ['DataMask','Water','NotBuilt','b2014','b2000','b1990','b1975'])
            curRes = processGHSLsummary(outFile, False)
            row[1:] = curRes
            cursor.updateRow(row)
    for f in allOutFiles:
        os.remove(f)
            
###Open GHSL file - 
# 1. Check if the file has had its projection defined
# 2. If the projection isn't Web Mercator, define it (NOT a re-projection) as web mercator
def openGHSL(gFile, ret = True):
    gPrj = arcpy.Describe(gFile).spatialReference
    if gPrj.name == "WGS_1984_Pseudo_Mercator":
        arcpy.DefineProjection_management (gFile, arcpy.SpatialReference(r"S:\GLOBAL\Global_Human_Settlement\BETA\index2.prj"))
    else:
        fubar=0
    if ret:
        return(Raster(gFile))
    return(1)

def getListVals(xx, start):
    return(sum([xx[x] for x in range(start,len(xx),5)]))

#Summarize the datamask by the builtup layer - For each category in the input builtup (3-6) tabulate the datamask values
#inMask = r"S:\GLOBAL\Global_Human_Settlement\BETA\FULL\DATAMASK\12\2880\2576_datamask.tif"
#inBuiltup = r"S:\GLOBAL\Global_Human_Settlement\BETA\FULL\MT\12\2880\2576_mt.tif"
def summarizeMaskByBuilt(inBuiltup, inMask):
    inB = arcpy.RasterToNumPyArray(inBuiltup)
    inM = arcpy.RasterToNumPyArray(inMask)
    outDict = {}
    for bType in [1,2,3,4,5,6]:
        #Create a mask for the current built up area
        curB = np.ma.getmask(np.ma.masked_equal(inB, bType))
        maskedM = curB * inM   #Multiply the data mask by the builtup mask
        #Tabulate the masked data
        y = np.bincount(np.concatenate(maskedM))
        #ii = np.nonzero(y)[0]    
        #xx = zip(ii,y[ii])
        for idx in range(1, 17):
            try:
                curVal = y[idx]
            except:
                curVal = 0
            if not idx in outDict.keys():
                outDict[idx] = [curVal]
            else:
                outDict[idx].append(curVal)
    return outDict

def selectGHSLtiles(shapeFL, tileOutlines, locationField, ghslProcess=True):
    outlinesFL = "GHSLOut"
    #1. Select all tiles in the ghslOutlines that intersect the input shape
    arcpy.MakeFeatureLayer_management(tileOutlines, outlinesFL)
    arcpy.SelectLayerByLocation_management(outlinesFL, 'INTERSECT', shapeFL)
    curTitles = []
    with arcpy.da.SearchCursor(outlinesFL,[locationField]) as cursor:
        for feat in cursor:
            curTitle = feat[0]   
            if ghslProcess:
                curTileSpl = curTitle.split("/")        
                curTitles.append("%s\\%s" % (curTileSpl[len(curTileSpl)-2], curTileSpl[len(curTileSpl)-1].replace("band_1_href","mt")))
            else:
                curTitles.append(curTitle.replace("Settlement\\BETA", "Settlement_Layer\\BETA"))               
    return curTitles

def processGHSLsummary(inFile, includeQuality=True):
    inD = np.genfromtxt(inFile, delimiter=',')    
    #Summarized number of pixels in each of the time periods 2014, 2000, 1996, 1975
    dataQuality = [inD[14:,3:7].sum(), inD[10:14,3:7].sum(), inD[6:10,3:7].sum(), inD[1:6,3:7].sum()]
    #Summarized built in 2014 only and built base (everything from 2000 backwards in time)
    builtPeru = [inD[1:,3].sum(), inD[1:,4:7].sum()]
    #Summarized epochs
    builtTotal = [inD[1:,3].sum(), inD[1:,4].sum(), inD[1:,5].sum(), inD[1:,6].sum()]    
    #Summarized Water and Not Built
    notBuilt = [inD[1:,2].sum(), inD[1:,1].sum()]
    results = notBuilt + builtTotal
    if includeQuality:
        results = results + dataQuality
    return(np.nan_to_num(results))
    
def addGUF(gufMap, curTiles, urbanLayer, curMap, gufFolder, gufSymbology, extentsSymbology):
    outlinesFL = "GHSLOut2"
    shapeFL = "InShape"
    
    #Open map and add urban layer
    mxd = arcpy.mapping.MapDocument(gufMap)
    mxd.saveACopy(curMap)
    mxd = arcpy.mapping.MapDocument(curMap)
    df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
    urbLyr = arcpy.mapping.Layer(urbanLayer)
    urbLyr.definitionQuery = 'NOT  "ExtentName" = \'0\''
    arcpy.ApplySymbologyFromLayer_management(urbLyr, extentsSymbology)
    arcpy.mapping.AddLayer(df, urbLyr)
    
    lyrCnt = 0
    for c in curTiles:
        c = os.path.join(gufFolder, c)
        outLayer = "GHSL_%s" % lyrCnt        
        #if not outLayer in [l.name for l in arcpy.mapping.ListLayers(mxd)]:
        result = arcpy.MakeRasterLayer_management(c, outLayer)
        layer = result.getOutput(0)
        arcpy.ApplySymbologyFromLayer_management(layer, gufSymbology)
        arcpy.mapping.AddLayer(df, layer)
        lyrCnt = lyrCnt +1
    df.extent = arcpy.Describe(urbanLayer).extent
    mxd.save()                 
    
  

#CODE to add GHSL raster to ArcMap
def addGHSL(urbanLayer, curMap, ghslFolder, datamaskFolder, lyrSymbology, extentsSymbology, ghslOutline):    
    inputGHSLMap = r"S:\GLOBAL\Global_Human_Settlement_Layer\GHSL_Footprint_Maps.mxd"
    outlinesFL = "GHSLOut2"
    shapeFL = "InShape"
    arcpy.MakeFeatureLayer_management(ghslOutline, outlinesFL)

    #Open map and add urban layer
    mxd = arcpy.mapping.MapDocument(inputGHSLMap)
    mxd.saveACopy(curMap)
    mxd = arcpy.mapping.MapDocument(curMap)
    df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
    urbLyr = arcpy.mapping.Layer(urbanLayer)
    urbLyr.definitionQuery = 'NOT  "ExtentName" = \'0\''
    arcpy.ApplySymbologyFromLayer_management(urbLyr, extentsSymbology)
    arcpy.mapping.AddLayer(df, urbLyr)
    
    #Find intersecting
    curTitles = selectGHSLtiles(urbanLayer, ghslOutline, 'Path', False)
    
    lyrCnt = 0

    for c in curTitles:
        c = os.path.join(ghslFolder, c)
        #Add MT layer
        outLayer = "%s_%s" % (os.path.basename(c), "mt")        
        #if not outLayer in [l.name for l in arcpy.mapping.ListLayers(mxd)]:
        result = arcpy.MakeRasterLayer_management(c, outLayer)
        layer = result.getOutput(0)
        arcpy.ApplySymbologyFromLayer_management(layer, lyrSymbology)
        arcpy.mapping.AddLayer(df, layer)

    df.extent = arcpy.Describe(urbanLayer).extent
    mxd.save()                 
    
def createGHSLexcel (ghslShp, outGHSL):
    #Write ghsl shape to workbook
    writeShapefileXLS(outGHSL, ghslShp, "GHSL", 'pop', 'D')
    workbook = openpyxl.load_workbook(outGHSL)
    dataSheet = workbook["GHSL"]
    #Add total built up columns to sheets
    dataSheet.cell(row=1, column=11).value = "built1990"
    dataSheet.cell(row=1, column=12).value = "built2000"
    dataSheet.cell(row=1, column=13).value = "built2015"
    for rowIdx in range(2, 20):
        dataSheet.cell(row=rowIdx, column=11).value = "=sum(I%s:J%s)" % (rowIdx, rowIdx)
        dataSheet.cell(row=rowIdx, column=12).value = "=sum(H%s:J%s)" % (rowIdx, rowIdx)
        dataSheet.cell(row=rowIdx, column=13).value = "=sum(G%s:J%s)" % (rowIdx, rowIdx)
    
    chartPage = workbook.create_sheet("Charts")
    extentChart = openpyxl.chart.BarChart()
    extentChart.type = "col"
    extentChart.style = 10
    extentChart.title = 'GHSL'
    extentChart.x_axis.title = "10 largest cities"
    extentChart.y_axis.title = 'Total Built Area'
    extentChart.width = 20
    extentChart.height = 15
    
    data = openpyxl.chart.Reference(dataSheet, min_col=10, max_col=13, min_row=1, max_row=11)
    cats = openpyxl.chart.Reference(dataSheet, min_col=3, min_row=2, max_row=11)
    extentChart.add_data(data, titles_from_data=True)
    extentChart.set_categories(cats)    
    chartPage.add_chart(extentChart, "A1")
    workbook.save(outGHSL)
    workbook.close()