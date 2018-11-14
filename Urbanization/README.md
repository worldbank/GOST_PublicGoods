# Urban Toolbox
This repository is the aggregation of a number of tools developed for and by the World Bank Group for classifying, quantifying and creating urban analytics. The GOSTRocks folder contains the library for esecuting these tasks, and the ExecutionScripts folder contains scripts for using that library.

## Frequently used tools
### European Union Urban Classification

``` python
import GOSTRocks.Urban.WorldPop_Library as wp
import GOSTRocks.arcpyMisc as arcpyMisc


lowDensVal = 300
highDensVal = 1500
lowDensPop = 5000
highDensPop = 50000
#These two can be the same if the cell size is 1 km2
tempPop = "C:/PATH/TO/POPULATIONRASTER.tif" 
tempDen = "C:/PATH/TO/POPULATIONDENSITYRASTER.tif"

wp.createUrbanClusters(tempPop, tempDen, lowDensVal, urbClstGrid)
wp.createUrbanClusters(tempPop, tempDen, highDensVal, hdClstGrid)
   
arcpy.env.workspace = popOutputFolder
wp.smoothClusters(tempDen, tempPop, hdClstGrid, hdClstGridSmooth, highDensPop, 16, False)
arcpy.CalculateStatistics_management(hdClstGridSmooth, "1", "1", "#", "OVERWRITE")

#Convert rasters to shapefiles
arcpyMisc.extractFootprints(urbClstGrid, 4999, urbShp)
arcpyMisc.extractFootprints(hdClstGridSmooth, 49999, hdShp)

straight = wp.summarizeWorldPop(urbClstGrid, tempPop, urbClstTbl, "ID")  
smoothed = wp.summarizeWorldPop(hdClstGridSmooth, tempPop, hdClstTbl, "ID")              
```
