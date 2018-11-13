
import sys, os, inspect
cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import GOSTRocks.arcpyMisc
import GOSTRocks.Urban.UrbanReviewv2
#import GOSTRocks.Urban.accessibilityMap
#import pandas as pd

#Calculate Puga index for each feature in a shapefile
import GOSTRocks.Urban.PugaIndex as puga
inShp = r"C:\temp\OECD\OECD_cities_sample2.shp"
puga.main(inShp, inShp.replace(".shp", "_PUGA.csv"), os.path.dirname(inShp))


'''
xx = GOSTRocks.arcpyMisc.getLargestPointsInPolygons(r"C:\Temp\TUN\GIS\ADMIN\TUN_2.shp", 
    r"S:\GLOBAL\POPULATION\Cities_D\GRUMP\global_settlement_points_v1_01.shp",
    "pop", "WB_ADM2_NA", "Schnm")
print xx
    
d = pd.DataFrame({'id' : pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
     'Long' :pd.Series([-77.036815, -77.019547, -77.036815, -77.019547], index=['a', 'b', 'c', 'd']),
     'Lat' : pd.Series([38.9093104, 38.9093104, 38.9126372, 38.9126372], index=['a', 'b', 'c', 'd']),
     'Weight' : pd.Series([1000., 2000., 30000., 4000.], index=['a', 'b', 'c', 'd'])})

xx = GOSTRocks.Urban.accessibilityMap.accessMap(d)
print xx.calcODMatrix()

countries = ['VNM']
for c in countries:
    GOSTRocks.Urban.UrbanReviewv2.calculateUrban(c, r"C:\Temp", cAdmin=False, mapCities=False, iUrban=False)
for c in ['ALB']: #'BHR','DZA',FAILED: 'ISR','PSE',
	GOSTRocks.Urban.UrbanReviewv2.calculateUrban(c, r"C:\Temp",
        finalizeProduct=False, mapCities=False, 
        cPopulation=True, cNTL=True, cGHSL=True, cGUF=True, cAI=False, cAdmin=False, 
        cPuga=True, cCompactness=True, cDiscontiguity=True)
'''
