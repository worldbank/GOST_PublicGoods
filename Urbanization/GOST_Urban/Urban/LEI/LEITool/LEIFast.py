# -*- coding: utf-8 -*-

import arcpy
from arcpy import env

iWorkspace=arcpy.GetParameterAsText(0)
oriShp=arcpy.GetParameterAsText(1)
newShp=arcpy.GetParameterAsText(2)
dist=arcpy.GetParameterAsText(3)



def LEIFast(iWorkspace, oriShp, newShp, dist):
    try:
        
        arcpy.AddMessage('Step 1/2. Open data')
        print('Step 1/2. Open data')
        dist=float(dist)
        env.workspace=iWorkspace
        ##Check new layer has field "LEI", if not then creat
        newLayer="new_layer"
        arcpy.MakeFeatureLayer_management(newShp, newLayer)
        desc=arcpy.Describe(newLayer)
        ifieldInfo=desc.fieldInfo
        index=ifieldInfo.findfieldbyname("LEI")
        if index==-1:
            arcpy.AddField_management(newLayer,"LEI","DOUBLE")
        ##Check old layer has field "MLEI", if not then creat and set a defult value 100
        oriLayer="oriLayer"
        arcpy.MakeFeatureLayer_management(oriShp, oriLayer)
        descOri=arcpy.Describe(oriLayer)

        arcpy.AddMessage('Step 2/2. Caculate LEI')
        print('Step 2/2. Caculate LEI')
        ##Compute the LEI index
        GeometryRows=arcpy.UpdateCursor(newLayer)
        for newFeature in GeometryRows:
            iFDifference=newFeature.shape.buffer(dist).difference(newFeature.shape)
            arcpy.SelectLayerByLocation_management(oriLayer,"INTERSECT",iFDifference)
            GeometryRows2=arcpy.SearchCursor(oriLayer)
            inAreaLEI=0
            inAreaMLEI=0
            error = 0
            for oldFeature in GeometryRows2:
	        try:
                    insideArea=iFDifference.intersect(oldFeature.shape,4)
                    inAreaLEI+=insideArea.area
		except:
		    error = 1
		    break
	    if error == 0:
	        newFeature.LEI = inAreaLEI/iFDifference.area*100
	    else:
	        newFeature.LEI = 999
            GeometryRows.updateRow(newFeature)
            del GeometryRows2
        del GeometryRows

        arcpy.AddMessage('Finished!')
        print('Finished!')
    
    except Exception as e:
        arcpy.AddMessage("Error: " + str(e.message))

		
		
if __name__=="__main__":

    """
    iWorkspace="F:\Projects\Liu\LEI\GuangZhou"
    oriShp="gz2005.shp"
    newShp="gz2010.shp"
    dist=200
    """

    LEIFast(iWorkspace, oriShp, newShp, dist)
