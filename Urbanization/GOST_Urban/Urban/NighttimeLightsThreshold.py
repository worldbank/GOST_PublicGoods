###############################################################################
# Extract ESA Landcover poitns NTL values
# Benjamin P Stewart, Apr 2015
###############################################################################

import time, sys, os, csv, re, random
from GOSTRocks.misc import *
from GOSTRocks.arcpyMisc import *
from GOSTRocks.listStuff import *

import arcpy
from arcpy.sa import *

arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

#Read in Countries from output threshold CSV
def attributeESApoints(cntry, outFeatures, outDB, ntlFile, lcRaster, ghslValu, globalOutlines):    
    tempLC = "%s/%s" % (outDB,"ESA2009")
    tempPts = "%s/%s" % (outDB,"tPts")
    tempAdmin = "%s/%s" % (outDB, "Admin1")
    
    arcpy.Select_analysis(in_features=globalOutlines, out_feature_class=tempAdmin, where_clause=""""ISO3" = '%s'""" % cntry)
    tPrint("***Created Admin File")
    
    arcpy.Clip_management(in_raster=lcRaster, rectangle="", out_raster=tempLC, in_template_dataset=tempAdmin, nodata_value="0", clipping_geometry="ClippingGeometry")
    tPrint("***Clipped ESA Globcover")
    
    arcpy.RasterToPoint_conversion(in_raster=tempLC, out_point_features=tempPts, raster_field="Value")
    renameField(tempPts, "grid_code", "Globcover")
    tPrint("***Converted to points")
    
    ExtractValuesToPoints(tempPts, ntlFile, outFeatures, "NONE", "VALUE_ONLY")
    renameField(outFeatures, "RASTERVALU", "NTL")
    tPrint("***Extracted NTL Values")
    fOut = outFeatures.replace("Globcover", "Globcover_GHSL")
    
    arcpy.sa.ExtractValuesToPoints(outFeatures, ghslValu, fOut)
    renameField(fOut, "RASTERVALU", "GHSL")
    arcpy.Delete_management(outFeatures)
    arcpy.Rename_management(fOut, outFeatures)
    

'''
Create an iterable range object that allows decimals as steps
'''
def getHistIndex(hIdx, val, maxVal=2000):
    for h in range(0, len(hIdx)):
        curH = hIdx[h]
        if curH > val:
            return(lastH)        
        lastH = h
    return(len(hIdx) -1)
        
def getHistPer(inD):
    tSum = listSum(inD)
    for hIdx in range(0,len(inD)):
        inD[hIdx] = inD[hIdx] / tSum
    return(inD)

def listCumSum(inD, rev=False):
    outD = [0] * len(inD)
    if not rev:
        outD[0] = inD[0]
        for idx in range(0,len(inD)):
            outD[idx] = outD[(idx - 1)] + inD[idx]
    else:
        outD[len(outD)-1] = inD[len(outD)-1]
        for idx in range(len(inD)-2, 0, -1):
            outD[idx] = outD[(idx + 1)] + inD[idx]
    return(outD)

def listAdd(d1, d2):
    outD = [0] * len(d1)
    for idx in range(0, len(d1)):
        outD[idx] = d1[idx] + d2[idx]
    return(outD)
'''
Calculate the bins of a 0.5 step histogram from the NTL rastervalu column
calcField="RASTERVALU"
urbField="grid_code"
urbDef=190
ntlMax = 1000
inFC = r"C:\Work\NTL\CityLights\ThreshTesting\SampledPoints.gdb\BLR_Globcover_RC2010"

'''
def calculateHistogram(inFC, ntlMax, calcField="NTL", urbField="NTL", ghsl=False):
    urbDef=190
    hDef = [x for x in drange(0,ntlMax,0.5)]
    urbVal = [0] * len(hDef)    #Histogram of urban
    nonVal = [0] * len(hDef)    #Histogram of non-urban    
    #create a search cursor for the current fc, grab the rasterValu, and 
    fcCur = arcpy.da.SearchCursor(inFC, [calcField,urbField], "%s > 0" % calcField)
    rowCnt = 0
    for row in fcCur:
        urban = row[1] == urbDef
        if ghsl:
            urban = row[1] > 127
        histIndex = getHistIndex(hDef, row[0])
        if urban:
            urbVal[histIndex] += 1
        else:
            nonVal[histIndex] += 1
        rowCnt += 1
        if rowCnt % 1000000 == 0:
            tPrint("Processed Row %s" % rowCnt)            

    #Convert the histograms to percentages
    perUrb = getHistPer(urbVal)
    perNonUrb = getHistPer(nonVal)

    #Calculate the cumulative sum for each of urban and non-urban for each threshold
    cumUrb = listCumSum(perUrb, rev=True)
    cumNonUrb = listCumSum(perNonUrb)
    tAcc = listAdd(cumUrb, cumNonUrb)
    threshIdx = listMax(tAcc, idx=True)
    threshold = hDef[threshIdx]

    return({'histDef':hDef, 'urbanHistogram':perUrb, 'nonUrbanHistogram':perNonUrb, 'threshold':threshold})

def calculateThreshold(curCountry, outCSV="S:/GLOBAL/Projects/CityLights/Data/Tabular/ThreshComp.csv", 
    outDB=r"S:\GLOBAL\Projects\CityLights\Data\Vector\SampledPoints.gdb", 
    ntlFile="S:/GLOBAL/NightLights/rad_cal/F2010_0111-20101209_rad_v4.avg_vis_Corrected.tif",
    histFolder = "S:/GLOBAL/Projects/CityLights/histograms",
    lcRaster = r"S:\GLOBAL\Landcover_D\Globcover_D\GLOBCOVER_L4_200901_200912_V2.3.tif",
    ghslValu = r"S:\GLOBAL\Global_Human_Settlement_Layer\BETA\AGG\T2014_300m_Mollweide.tif",
    globalOutlines = "S:/GLOBAL/ADMIN/WB-2014_subadmin_D/Polygons/Admin0_Polys.shp"):
    #Get a list of countries with thresholds calculated
    processedCountries = []
    with open(outCSV, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            processedCountries.append(row[0])
            if row[0] == curCountry:
                tPrint("***Already has a calculated threshold")
                return(row[1])
    processedCountries = processedCountries[1:]  
    #If country is not processed the calculate the threshold
    if not curCountry in processedCountries:
        threshWriter = open(outCSV, 'a')
        outFeatures = "%s/%s_Globcover_RC2010" % (outDB,curCountry)
        finalFeatures = "%s/%s_Globcover_RC2010_Sampled.shp" % (outDB, curCountry)
        if not arcpy.Exists(outFeatures):
            #Create output ESA points shapefile with NTL attributed			
            attributeESApoints(curCountry, outFeatures, outDB, ntlFile, lcRaster, ghslValu, globalOutlines)
            tPrint("***Extracted RadCal 2010 NTL data")
        else:
            tPrint("***%s already has attributed ESA points" % curCountry)
        xx = calculateHistogram(outFeatures, 2000)    
        writer = open("%s/%s_hist.csv" %(histFolder, curCountry), 'w')
        writer.write("NTLVal,Urban,NonUrban\n")
        for xIdx in range(0,len(xx['histDef'])):
            writer.write('%s,%s,%s\n' %(xx['histDef'][xIdx],xx['urbanHistogram'][xIdx],xx['nonUrbanHistogram'][xIdx]))
        writer.close()
        tPrint("***%s processed with a threshold of %s" % (curCountry, xx['threshold']))
        threshWriter.write("%s,%s\n" % (curCountry, xx['threshold']))
        threshWriter.close() 
        return(xx['threshold'])
    else:
        tPrint("***Already has a calculated threshold")    

def summarizeNTLUrbanAreas(curCountry, inShape, id, sumVars=["SUM"]):    
    allFiles = []
    inTitles = [id]
    inFiles = getNTLFiles()
    #Calculate threshold for current country using both raw and rc 
    threshCSV="S:/GLOBAL/Projects/CityLights/Data/Tabular/ThreshCompRaw.csv"
    threshDB=r"S:\GLOBAL\Projects\CityLights\Data\Vector\SampledPoints_Raw.gdb"
    ntl=r"S:\GLOBAL\NightLights\2009\F162009.v4b_web.stable_lights.avg_vis_ElvidgeCorrected_gt3.tif"
    rawThresh = calculateThreshold(curCountry, outCSV = threshCSV, outDB=threshDB, ntlFile=ntl)
    
    threshCSV="S:/GLOBAL/Projects/CityLights/Data/Tabular/ThreshComp.csv"
    threshDB=r"S:\GLOBAL\Projects\CityLights\Data\Vector\SampledPoints.gdb"
    ntl=r"S:\GLOBAL\NightLights\rad_cal\F2010_0111-20101209_rad_v4.avg_vis_Corrected.tif"
    rcThresh = calculateThreshold(curCountry, outCSV = threshCSV, outDB=threshDB, ntlFile=ntl)
    fCount = 0
    #Summarize all the NTL files for the current inShape
    for f in inFiles:
        #Create titles for the ntl files for titling output columns
        fName = os.path.basename(f)[:7]
        if "_" in fName:
            fName = fName.split("_")[0]
            fName = fName.replace("F", "RC")       
        inTitles.append("%s_URBAN" % fName)
        inTitles.append("%s_TOTAL" % fName)
        
        #Determine the appropriate threshold based on the filename
        ntlType = os.path.basename(f)[7:8]
        if ntlType == ".":
            curThresh = rawThresh
        else:
            curThresh = rcThresh
        
        tPrint("Processing %s for %s" % (curThresh, f))
        #Run the analysis on both urban and rural lights
        outDBFUrban = "C:/Temp/%s_%s_Urb.dbf" % (curCountry, fCount)
        outDBFRural = "C:/Temp/%s_%s_Rur.dbf" % (curCountry, fCount)
        allFiles.append(outDBFUrban)
        allFiles.append(outDBFRural)
        if not os.path.exists(outDBFUrban):
            summarizeSum(f, curThresh, inShape, id, outDBFUrban)        
            summarizeSum(f, curThresh, inShape, id, outDBFRural, False)        
        fCount += 1
    
    if id in ['FID', 'OID']:
        id = "%s_" % id
        
    res = summarizeDbf(allFiles, id, -1)
    for f in allFiles:
        os.remove(f)
    return({"Results":res, "Titles":inTitles})

    
def summarizeArea(inRaster, thresh, inShape, id, outDBF):
    inR = Raster(inRaster)
    tRaster = inR > float(thresh)
    ZonalStatisticsAsTable(inShape, id, tRaster, outDBF)

def summarizeSum(inRaster, thresh, inShape, id, outDBF, gt=True):
    inR = Raster(inRaster)
    if gt:
        tRaster = (inR > float(thresh)) * inR
    else:
        tRaster = (inR < float(thresh)) * inR
    
    ZonalStatisticsAsTable(inShape, id, tRaster, outDBF)
    