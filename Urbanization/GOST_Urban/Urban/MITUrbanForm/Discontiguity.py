# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Discontiguity
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import arcpy, sys, os
sys.path.append(r"C:\Users\wb411133\Box Sync\AAA_BPS\Code\GOST")
from GOSTRocks.misc import tPrint
from GOSTRocks.arcpyMisc import tryAddField


def calcDiscontiguityFromGUF(extents, gufRaster, tempFolder, fieldPref="GUF"):
    tryAddField(extents, "%sCont" % fieldPref, "FLOAT")
    
    with arcpy.da.UpdateCursor(extents, ["OID@", "%sCont" % fieldPref]) as cursor:
        for featRow in cursor:
            tRaster = os.path.join(tempFolder, "pop_%s" % featRow[0])
            tPoly = os.path.join(tempFolder, "pop_%s.shp" % featRow[0])
            tLayer = "pop_%s" % featRow[0]            
            
            arcpy.MakeFeatureLayer_management(extents, tLayer, '"FID" = %s' % featRow[0])
            #Use the current feature to extract the current raster
            try:
                arcpy.Clip_management(gufRaster, '#', tRaster, tLayer, '0', 'ClippingGeometry', 'MAINTAIN_EXTENT')
                arcpy.RasterToPolygon_conversion(tRaster, tPoly, "NO_SIMPLIFY")
                featRow[1] = calculateDiscontiguity(tPoly, "square_kilometer")
            except:
                print ("Something went wrong with feature %s" % featRow[0])
                featRow[1] = 0
            cursor.updateRow(featRow)              

def calculateDiscontiguity(built_up_area, unit):
    fields = arcpy.ListFields(built_up_area)
    fields = [i.name for i in fields]
    if "AREA" not in fields:
        arcpy.AddField_management (built_up_area, "AREA","DOUBLE")

    dic_unit = {'square_mile':'squaremile','square_meter':'squaremeter',\
               'square_foot':'squarefoot','square_miile':'squaremile',\
               'square_kilometer':'squarekilometer'}   
    arcpy.CalculateField_management (built_up_area, "AREA","!shape.area@{0}!".format(dic_unit[unit]),'PYTHON')    
    cursor=arcpy.da.SearchCursor(built_up_area,["AREA","FID"])
    total_area = 0
    raw_index = 0
    for row in cursor:        
        total_area = total_area + row[0]
        sum_Area_fCs_smaller_than_current_FC = 0
        inner_cursor = arcpy.da.SearchCursor(built_up_area,"AREA",'"AREA"<={0} AND "FID"<>{1}'.format(row[0],row[1]))
        for inner_row in inner_cursor:
            sum_Area_fCs_smaller_than_current_FC = sum_Area_fCs_smaller_than_current_FC + inner_row[0]
        del inner_cursor
        raw_index = raw_index + sum_Area_fCs_smaller_than_current_FC

    del cursor

    discontiguity=raw_index/total_area
    return(discontiguity)
    
    
    
    