################################################################################
# Launch Urban Analysis
# Benjamin P. Stewart, September 2018
# Purpose: Batch launch the municipal urban analysis
'''
import os
from shutil import copyfile

inFolder = r"Q:\WORKINGPROJECTS\CityScan\Data\CityExtents"

for root, dirs, files in os.walk(inFolder):
    if os.path.basename(root) == "Maps":
        for f in files:
            copyfile(os.path.join(root, f), os.path.join(r"C:\temp\Maps", f))
'''
################################################################################

import sys, os, inspect, logging
cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import GOSTRocks.Urban.UrbanAnalysis
import GOSTRocks.misc

logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.INFO)

'''
#Run landcover summaries for Balikpapan
inAOI = r"Q:\WORKINGPROJECTS\Indonesia_GBDx\BalikPapan_AOI.shp"
ua = GOSTRocks.Urban.UrbanAnalysis.urbanAnalysis(inShape = inAOI, urbanFolder = inAOI.replace(".shp", ""), gridSize=250)
inFile = r"Q:\WORKINGPROJECTS\Indonesia_GBDx\Balikpapan_GBDx\103001007E75D600\lulc\058350570010_01_assembly_clip_LULC.tif"
inFile2 = r"Q:\WORKINGPROJECTS\Indonesia_GBDx\Balikpapan_GBDx\104001003E70C400\lulc\058350569010_01_assembly_clip_LULC.tif"
ua.summarizeLULC(inFile , ua.lcFile.replace(".csv", "_GBDx_LULC_1.csv"))
ua.summarizeLULC(inFile2, ua.lcFile.replace(".csv", "_GBDx_LULC_2.csv"))
''' 

#CityScan Analysis
inFolder = r"Q:\WORKINGPROJECTS\CityScan\Data\CityExtents"
baseFloodFolder = r"Q:\GLOBAL\HYDRO\SSBN_Flooding"
'''
        ["Conotou","benin"],
        ["Banjul","gambia"],
        ["Douala","cameroon"],
        ["Freetown","sierra_leone"],
        ["Kampala","uganda"],
        ["Kigali","rwanda"],
        ["Monrovia","liberia"],
        ["Nairobi","kenya"],
        ["Zanzibar","tanzania"],
        ["Addis","ethiopia"],
        ["Bamako","mali"],
        ["Kinshasa","democratic_republic_congo"],cd 
        ["Abidjan","cote_d_ivoire"],
'''
data = [
        ["Dar","tanzania"]
       ]

urbanParams = GOSTRocks.misc.getUrbanParams()
cCnt = 0
for curData in data:
    cCnt += 1
    c = curData[0]
    floodFolder = os.path.join(baseFloodFolder, curData[1])
    curShp = os.path.join(inFolder, "%s_AOI.shp" % c)
    if os.path.exists(curShp):
        logging.info("Processing %s" % curShp)
        outFolder = curShp.replace(".shp", "")
        if not os.path.exists(outFolder):
            os.mkdir(outFolder)
        ua = GOSTRocks.Urban.UrbanAnalysis.urbanAnalysis(inShape=curShp, urbanFolder=outFolder, gridSize=250)               
        ua.summarizeElevation()        
        ua.summarizeOSM()
        ua.prepMappingData()
        '''
        print(ua.summarizeData())
        ua.clipGHSL()
        ua.MarketAccess()
        ua.summarizeLandcover(urbanParams['afriCover20'])
        ua.summarizeGHSL()
        ua.searchForImagery()             
        ua.summarizeNDVI()        
        ua.summarizeBAI()
        ua.summarizeSSBN(floodFolder)
        ###This is the shitty stuff to use arcpy to create mapping documents
        ua = GOSTRocks.Urban.UrbanAnalysis.urbanAnalysis_arcpy(inShape=curShp, urbanFolder=outFolder, gridSize=250) 
        #ua.prepareMaps()
        ua.generateMaps()        
        '''
    else:
        logging.info("No extent present for %s" % curShp)