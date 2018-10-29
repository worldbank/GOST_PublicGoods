# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Creating Urban and High Density clusters from gridded population data
# By : Olivier Draily, DG REGIO, 2015-06-25
# Modified: Benjamin P Stewart, Sept, 2015
# Usage: Modify the local variables in this script and then run it (ArcGIS Desktop and Python must be installed on the computer).
#   Purpose is to create high density urban clusters and urban cluster above minimum
#   density and total population
#
# This script produces :
# - A raster of Urban Clusters (URB_CLST_GR) with the total population of each cluster.
# - A raster of High Density Clusters (HDENS_CLST_GR) with the total population of each cluster. 
# - Clusters are contiguous cells with a minimum density (ARGUMENT) and population (ARGUMENT).
# ---------------------------------------------------------------------------

from arcpy.sa import *
import arcpy, string, os, sys, traceback, time, csv
import numpy as np

def mapPop(inMap, outMap, hdGrid, urbGrid, hdStyle, urbStyle, extents, extentsSymbology):
    mxd = arcpy.mapping.MapDocument(inMap)
    mxd.saveACopy(outMap)
    mxd = arcpy.mapping.MapDocument(outMap)    
    df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
    
    #Add the urban grid to the map
    result = arcpy.MakeRasterLayer_management(urbGrid)
    urbLyr = result.getOutput(0)
    arcpy.ApplySymbologyFromLayer_management(urbLyr, urbStyle)
    arcpy.mapping.AddLayer(df, urbLyr)
    
    #Add the high density urban grid to the map
    resultHD = arcpy.MakeRasterLayer_management(hdGrid)
    hdLyr = resultHD.getOutput(0)
    arcpy.ApplySymbologyFromLayer_management(hdLyr, hdStyle)
    arcpy.mapping.AddLayer(df, hdLyr)
    
    #Add the extents to the map
    extentLyr = arcpy.mapping.Layer(extents)
    extentLyr.definitionQuery = 'NOT  "ExtentName" = \'0\''
    arcpy.ApplySymbologyFromLayer_management(extentLyr, extentsSymbology)
    arcpy.mapping.AddLayer(df, extentLyr)
    
    mxd.save()
    
    

##Runs Zonal statistics on the Urban Grid - each unique number in the Urban Grid represents
#   a unique feature (the total population of the feature)
#<RETURNS> - an array of 
#   1. The populations of the individual features in the urban grid
#   2. The total population
#   3. The population of the largest single feature
def summarizeWorldPop(urbGrid, popGrid, outTable, idVal="VALUE"):
    ZonalStatisticsAsTable(Int(urbGrid), idVal, popGrid, outTable)    
    allValues = []
    for r in arcpy.da.SearchCursor(outTable, ['SUM']):        
        allValues.append(r[0])       
    allVals = np.array(allValues)     
    return([allValues, np.sum(allVals), np.max(allVals)])

##Creates an urban cluster based on the population density grid and the population grid (these two grids are the same with the WorldPop data)
def createUrbanClusters(popDensGrid, popGrid, densVal, outGrid):
    #Classify below threshold as NoData, and above threshold as 1
    clustGrid = Reclassify(popDensGrid, "Value", "0 %s NODATA;%s 1000000 1" % (densVal, densVal), "DATA") 
    #Cluster contiguous features together
    groupedGrid = RegionGroup(clustGrid, "EIGHT", "WITHIN", "ADD_LINK", "")
    #Convert unique values grouped raster into a grid where unique values are the feature's total population
    sumUrbCluster = ZonalStatistics(groupedGrid, "VALUE", popGrid, "SUM", "DATA")
    #Convert 0 to threshold as NoData 
    finalGrid = Reclassify(sumUrbCluster, "VALUE", "0 %s NODATA" % densVal, "DATA")
    finalGrid.save(outGrid)    
    return(1)

#Cluster smoothing - This is applied to the High Density Clusters
def smoothClusters(popDensGrid, popGrid, clusterGrid, outGrid, highDensPop, minDiff=16, saveInterim=False):
    # Creation of R0_MASK0 where all values (population >= 0) are converted to 0. It will be used later in the smoothing process.
    r0Test = Raster(popGrid) > 0
    r0Mask = Reclassify(r0Test, "VALUE", "0 NoData; 1 0", 'DATA')
    # Creation of a raster where gaps are filled
    clusterMajority = MajorityFilter(clusterGrid, "FOUR", "MAJORITY")
    # Majority Filter filled gaps inside clusters, but it also deleted a lot of small clusters. We must recover them.
    clusterMajority0 = Reclassify(clusterMajority, "VALUE", "NoData 0; 1 10000000 1", "DATA") # NODATA -> 0 to make sa.Plus
    clusterFinal = Plus(clusterMajority0, clusterGrid)
    # Creation of a raster where only the clusters (now filled) with a minimum population of 50 000 inhabitants are kept
    # The extent of R7_PLUS is correct for this, but high-density cells are no more clustered => RegionGroup_sa again
    # => a new mask is needed to avoid assigning a number outside the high-density cells
    clustMask = Reclassify(clusterFinal, "Value", "0 NODATA;1 10000000 1", "DATA") # mask creation
    arcpy.env.mask = clustMask
    clustMaskGrouped = RegionGroup(clustMask, "FOUR", "WITHIN", "NO_LINK", "")
    arcpy.env.mask = ""
    clustPopZonal = ZonalStatistics(clustMaskGrouped, "VALUE", popGrid, "SUM", "DATA")
    if saveInterim:
        clustPopZonal.save("clustPopZonal_L48.tif")
    clustPopGTPop = Reclassify(clustPopZonal, "VALUE", "0 %s NODATA;%s 50000000 1" % (highDensPop, highDensPop), "DATA")
    if saveInterim:
        clustPopGTPop.save("clustPopGTPop_L53.tif")
    # To smooth the clusters, they must have a unique number => regiongroup
    clustPopUnq = RegionGroup(clustPopGTPop, "FOUR", "WITHIN", "ADD_LINK", "")
    # To smooth the clusters, the background must also be 0
    clustPopGTPop0 = Reclassify(clustPopGTPop, "VALUE", "NODATA 0; 1 1", "DATA")
    clustPopGTPopCon = Con(clustPopGTPop0, "0", clustPopUnq, "VALUE=0") # combine 0 from R13_GT50K_RECL and the unique ID from R12_GT50K_RGG3
    # Small bays in the high-density clusters are smoothed using the majority rule iteratively.
    # The majority rule means that if at least five out of the eight cells surrounding a cell belong to the same high-density cluster, it will be added.
    sumDiff = 1200
    i = 1
    R15_HDC_SMTH = clustPopGTPopCon
    if saveInterim:
        R15_HDC_SMTH.save("R15_HDC_SMTH_L66.tif")
    while (abs(sumDiff) > minDiff and i < 50):
        R16_Majority = MajorityFilter(R15_HDC_SMTH, "EIGHT", "MAJORITY")  # "EIGHT"
        R17_RECL_MAJ = Reclassify(R16_Majority, "VALUE", "0 0;1 1000000 1", "DATA")
        R18_RECL_HDC = Reclassify(R15_HDC_SMTH, "VALUE", "0 0;1 1000000 1", "DATA")
        R19_PLUS = Plus(R17_RECL_MAJ, R18_RECL_HDC)
        R20_PLUS_RECL = Reclassify(R19_PLUS, "VALUE", "0 NODATA;1 2 1", "DATA")
        R21_RGG = RegionGroup(R20_PLUS_RECL, "FOUR", "WITHIN", "NO_LINK", "")
        R22_IsNull = IsNull(R21_RGG)
        R23_CON = Con(R22_IsNull, "0", R21_RGG, "VALUE =1")
        R23_SMOOTHING = Plus(r0Mask, R23_CON)
        R24_MINUS = Minus(R23_SMOOTHING, R15_HDC_SMTH)
        ZonalStatisticsAsTable(r0Mask, "Value", R24_MINUS, "DIFF_R_S_TEMP", "DATA", "SUM")
        cur = arcpy.da.SearchCursor("DIFF_R_S_TEMP", ['SUM'])
        for row in cur:
            sumDiff = row[0]
        del cur, row
        #if i % 5 == 0:
        if saveInterim:
            R23_SMOOTHING.save("%s%s" % (outGrid, i))
        R15_HDC_SMTH = R23_SMOOTHING
        i = i + 1
        print i
        
        
    M_Clc_Ppl_GR = Reclassify(R15_HDC_SMTH, "Value", "0 NODATA;1 90000 1", "DATA")
    REG_GP_POPL = RegionGroup(M_Clc_Ppl_GR, "EIGHT", "WITHIN", "ADD_LINK", "")
    if saveInterim:
        REG_GP_POPL.save("REG_GP_POPL")
    zonalPop = ZonalStatistics(REG_GP_POPL, "VALUE", popGrid, "SUM", "DATA")
    zonalPop = Int(zonalPop)
    zonalPop.save(outGrid)


    