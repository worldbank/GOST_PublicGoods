###############################################################################
# Calculate Puga Index from population data
# Benjamin P. Stewart and __author__ = 'SPIJKERM'
# Purpose: Not sure ... will fill up later
###############################################################################
__author__ = 'GOST and SPIJKERM'

import os, sys, inspect
cmd_folder = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
from GOSTRocks.misc import tPrint
from GOSTRocks.arcpyMisc import summarizeDbf
from GOSTRocks.arcpyMisc import writeDict
import arcpy
from arcpy import env
from arcpy.sa import *

#Environment settings
arcpy.env.overwriteOutput = True
arcpy.env.projectCompare  = "Full"
arcpy.CheckOutExtension("spatial")

def main(cityExtent, outputFile, tempFolder,
        pop75 = r"Q:\GLOBAL\POP&DEMO\GHS\PopData\GHS-POP-LDS1975-GLOBE-R2015A_1975_54009_1000.tif",
        pop90 = r"Q:\GLOBAL\POP&DEMO\GHS\PopData\GHS-POP-LDS1990-GLOBE-R2015A_1990_54009_1000.tif",
        pop00 = r"Q:\GLOBAL\POP&DEMO\GHS\PopData\GHS-POP-LDS2000-GLOBE-R2015A_2000_54009_1000.tif",
        pop15 = r"Q:\GLOBAL\POP&DEMO\GHS\PopData\GHS-POP-LDS2014-GLOBE-R2015A_2015_54009_1000.tif"        
    ):
    #Local variables
    kernelFiles = os.path.join(os.path.dirname(__file__), "kernel%s_%s.txt")
    
    pop75_10 = calculatePuga(pop75, cityExtent, kernelFiles, 10, "Pop75", temp=tempFolder)    
    pop75_5 = calculatePuga(pop75, cityExtent, kernelFiles, 5, "Pop75", temp=tempFolder)
    pop90_10 = calculatePuga(pop90, cityExtent, kernelFiles, 10, "pop90", temp=tempFolder)
    pop90_5 = calculatePuga(pop90, cityExtent, kernelFiles, 5, "pop90", temp=tempFolder)
    pop00_10 = calculatePuga(pop00, cityExtent, kernelFiles, 10, "pop00", temp=tempFolder)
    pop00_5 = calculatePuga(pop00, cityExtent, kernelFiles, 5, "pop00", temp=tempFolder)
    pop15_10 = calculatePuga(pop15, cityExtent, kernelFiles, 10, "pop15", temp=tempFolder)
    pop15_5 = calculatePuga(pop15, cityExtent, kernelFiles, 5, "pop15", temp=tempFolder)
    l=[pop75_10, pop75_5, pop90_10, pop90_5, pop00_10, pop00_5, pop15_10, pop15_5]
    dbfList = [item for sublist in l for item in sublist]
    dbfHeaders = ["FID"]
    for x in dbfList:
        dbfHeaders.append("SUM_" + os.path.basename(x))
        dbfHeaders.append("MEAN_" + os.path.basename(x))    
        dbfHeaders.append("STD_" + os.path.basename(x))    
    
    allDbf = summarizeDbf(dbfList, "FID_", [['SUM','MEAN','STD']] * len(dbfList))    
    writeDict(allDbf, outputFile, dbfHeaders)
    #for f in dbfList:
    #    arcpy.Delete_management(f)
    
    
def calculatePuga(popFile, cityExtent, kernelFiles, distance, filePref="Pop", temp="C:/Temp", verbose=False):
    for f in [temp]:
        if not os.path.exists(f):
            os.makedirs(f)
        
    #Temporary files
    cityLS = os.path.join(temp, "%s_Rast.tif" % filePref)
    totalPop = os.path.join(temp, "%s_totalPop.tif" % filePref)

    focalPopulation = os.path.join(temp, "%s_focalPop_%s.tif" % (filePref, distance))
    peopleKernel1 = os.path.join(temp, "%s_People%skernel_1.tif" % (filePref, distance))
    peopleKernel2 = os.path.join(temp, "%s_People%skernel_2.tif" % (filePref, distance))
    pugaDistance = os.path.join(temp, "%s_Puga_%s.tif" % (filePref, distance))
    pugaKernel1 = os.path.join(temp, "%s_PugaKernel1_%s.tif" % (filePref, distance))
    pugaKernel2 = os.path.join(temp, "%s_PugaKernel2_%s.tif" % (filePref, distance))
    pugaDistanceDbf = os.path.join(temp, "%s_Puga_%s.dbf" % (filePref, distance))
    pugaKernel1Dbf = os.path.join(temp, "%s_PugaKernel1_%s.dbf" % (filePref, distance))
    pugaKernel2Dbf = os.path.join(temp, "%s_PugaKernel2_%s.dbf" % (filePref, distance))
    toDelete = [focalPopulation, peopleKernel1, peopleKernel2, pugaDistance, pugaKernel1, pugaKernel2]
    if not arcpy.Exists(pugaDistanceDbf):        
        #Clip the Landscan data to the city extent            
        arcpy.Clip_management(popFile, "", cityLS, cityExtent, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
        # calculate total population within the city
        ZonalStatistics(cityExtent, "FID", cityLS, "SUM").save(totalPop)
        if verbose:
            tPrint("Ran Zonal Stats to get city population")

        #Calculate people living in 10KM radius of cell i without discount factor, with e^-0.5d and e^-1d
        FocalStatistics (cityLS, NbrAnnulus(0,distance, "CELL"), "SUM", "").save(focalPopulation)    
        #Calculate people living in 10KM radius of cell i with distance decay functions - e^-0.5d and e^-1d
        FocalStatistics (cityLS, NbrWeight(kernelFiles % (distance, 1)), "SUM", "").save(peopleKernel1)
        FocalStatistics (cityLS, NbrWeight(kernelFiles % (distance, 2)), "SUM", "").save(peopleKernel2)
        if verbose:
            tPrint("Calculated Focal Statistics")

        (Raster(cityLS)/Raster(totalPop) * Raster(focalPopulation)).save(pugaDistance)
        (Raster(cityLS)/Raster(totalPop) * Raster(peopleKernel1)).save(pugaKernel1)
        (Raster(cityLS)/Raster(totalPop) * Raster(peopleKernel2)).save(pugaKernel2)
        if verbose:
            tPrint("Converted population to puga numbers")

        # Create Table of Puga outcomes
        ZonalStatisticsAsTable(cityExtent, "FID", pugaDistance, pugaDistanceDbf, "DATA", "ALL")
        ZonalStatisticsAsTable(cityExtent, "FID", pugaKernel1, pugaKernel1Dbf, "DATA", "ALL")
        ZonalStatisticsAsTable(cityExtent, "FID", pugaKernel2, pugaKernel2Dbf, "DATA", "ALL")
        toDelete.append(pugaDistance)
        toDelete.append(pugaKernel1)
        toDelete.append(pugaKernel2)
        if verbose:
            tPrint("Finished calculating Puga Index")
    for f in toDelete:
        arcpy.Delete_management(f)
    tPrint("Calculated Puga for %s and distance %s" % (filePref, distance))
    return([pugaDistanceDbf, pugaKernel1Dbf, pugaKernel2Dbf])

if __name__ == "__main__":
    main(cityExtents = sys.argv(1))