import sys, os, inspect, logging
import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon
import fiona
from rtree import index

#cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])))
cmd_folder = r"C:\Users\WB411133\OneDrive - WBG\AAA_BPS\Code\Code\Github\GOST"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import GOSTRocks.Urban.UrbanAdminComparison
import GOSTRocks.rasterMisc as rM
import GOSTRocks.misc as misc

urbanParams = misc.getUrbanParams()
adminFile = r"Q:\GLOBAL\ADMIN\GADM\v36_201805\gadm36ADM2.shp"
adminFile = r"Q:\GLOBAL\ADMIN\GADM\v36_201805\gadm36.shp"
inD = gpd.read_file(adminFile)
outFolder = r"Q:\GLOBAL\URBAN\GADM_Urban_Definition"
admCode = 'UID'

#Countries to process
inCountries = ["Afghanistan","Bangladesh","Burkina Faso","Chad","Ecuador","El Salvador","Fiji","Guinea","Haiti","Honduras","Madagascar","Mali","Mongolia","Myanmar","Mozambique","Niger","Pakistan","Papua New Guinea","Peru","Philippines","Senegal","Somalia","Thailand","Ukraine"]
for c in ['Somalia']:#inD.NAME_0.unique(): 
    curD = inD[inD.NAME_0 == c]        
    admFile = os.path.join(outFolder, "%s_adm2.shp" % c)
    euroStatFile = os.path.join(outFolder, "%s_GADM_lowestAdmin.csv" % c)
    ghslStatFile = os.path.join(outFolder, "%s_GADM_GHSL.csv" % c)
    
    if not os.path.exists(admFile):
        curD.to_file(admFile)        
    #Summarize GHSL
    ghslRes = rM.zonalStats(admFile, urbanParams['ghslVRT'], 
                bandNum=1, reProj = True, minVal = '',verbose=False , 
                rastType='C', unqVals=[0,1,2,3,4,5,6])
    ghslFinal = pd.DataFrame(ghslRes, columns = ['NoData','Water','NotBuilt','b00_14','b90_00','b75_90','bPre75'])
    curD_ghsl = pd.concat([curD, ghslFinal], axis=1)
    curD_ghsl.to_csv(ghslStatFile)
        
    #Run Eurostat urbanization methodology
    try:
        if not os.path.exists(euroStatFile):
            unionedFile = os.path.join(outFolder, "%s_unioned_areas_population.csv" % c)
            tabledOutput = unionedFile.replace(".csv", "_%s_Pivotted.csv" % c)
            newAdmin = adminFile.replace(".shp", "_noNull.shp")
            xx = GOSTRocks.Urban.UrbanAdminComparison.compareAreas(c, admFile, outFolder)
            unionedFile = xx.calculateVectorIntersection()
            tabulatedResults = xx.tabulateVectorResults(admCode=admCode)
            tabulatedResults.to_csv(euroStatFile)
    except:
        logging.warning("Could not process %s" % c)



