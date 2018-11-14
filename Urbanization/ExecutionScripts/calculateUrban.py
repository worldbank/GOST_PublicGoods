
import sys, os, inspect
cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import GOSTRocks.arcpyMisc
import GOSTRocks.Urban.UrbanReviewv2

#Calculate Puga index for each feature in a shapefile
import GOSTRocks.Urban.PugaIndex as puga
inShp = r"C:\temp\OECD\OECD_cities_sample2.shp"
puga.main(inShp, inShp.replace(".shp", "_PUGA.csv"), os.path.dirname(inShp))


'''
for c in ['ALB']: #'BHR','DZA',FAILED: 'ISR','PSE',
	GOSTRocks.Urban.UrbanReviewv2.calculateUrban(c, r"C:\Temp",
        finalizeProduct=False, mapCities=False, 
        cPopulation=True, cNTL=True, cGHSL=True, cGUF=True, cAI=False, cAdmin=False, 
        cPuga=True, cCompactness=True, cDiscontiguity=True)
'''
