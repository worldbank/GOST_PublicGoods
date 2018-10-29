# -*- coding: utf-8 -*-

import arcpy
from arcpy import env

iWorkspace=arcpy.GetParameterAsText(0)
oriShp=arcpy.GetParameterAsText(1)
newShp=arcpy.GetParameterAsText(2)
dist=arcpy.GetParameterAsText(3)



def LEI(iWorkspace, oriShp, newShp, dist):

    try:
        dist=float(dist)
        env.workspace=iWorkspace

        arcpy.AddMessage('Step 1/4. Open data')
        print('Step 1/4. Open data')

        #Get newLayer and check new layer has field "MLEI"&"LEI", if not then creat
        newLayer="new_layer"
        arcpy.MakeFeatureLayer_management(newShp, newLayer)
        desc=arcpy.Describe(newLayer)
        ifieldInfo=desc.fieldInfo
        index=ifieldInfo.findfieldbyname("LEI")
        if index==-1:
            arcpy.AddField_management(newLayer,"LEI","DOUBLE")
        #Get oriLayer/oldLayer
        oriLayer="oriLayer"
        arcpy.MakeFeatureLayer_management(oriShp, oriLayer)
        descOri=arcpy.Describe(oriLayer)

        #Create symmetrical difference layer to store the temple data
        if ".gdb" in iWorkspace:
            temple = "LEITemple"
        else:
            temple = "LEITemple.shp"
        templeLayer = "templeLayer"
        spatial_reference = arcpy.Describe(oriShp).spatialReference
        try:
            arcpy.CreateFeatureclass_management(iWorkspace, temple, "POLYGON", "", "DISABLED", "DISABLED", spatial_reference)
        except:
            arcpy.DeleteFeatures_management(temple)
        arcpy.MakeFeatureLayer_management(temple,templeLayer)
        if ".gdb" in iWorkspace:
            arcpy.AddField_management(templeLayer, "ID", "SHORT", "", "", "", "", "NULLABLE")

        arcpy.AddMessage('Step 2/4. Create difference between new polygon and its buffer')
        print('Step 2/4. Create difference between new polygon and its buffer')
     
        #Create the difference polygon between polygon and its buffer
        #Insert the polygon into temple layer
        tCursor = arcpy.InsertCursor(templeLayer)
        newCursor = arcpy.SearchCursor(newLayer)
        for newFeature in newCursor:
            iFDifference=newFeature.shape.buffer(dist).difference(newFeature.shape)
            row = tCursor.newRow()
            row.setValue("Shape",iFDifference)
            if ".gdb" in iWorkspace:
                row.setValue("ID",newFeature.OBJECTID)
            else:
                row.setValue("ID",newFeature.FID)
            tCursor.insertRow(row)
        del tCursor
        del newCursor
        del newFeature
        del row

        arcpy.AddMessage('Step 3/4. Caculate LEI')
        print('Step 3/4. Caculate LEI')
    

        #Compute the LEI index and creat a dictionary
        leiDict = dict()
        tCursor = arcpy.SearchCursor(templeLayer)
        for templeFeature in tCursor:
            arcpy.SelectLayerByLocation_management(oriLayer,"INTERSECT",templeFeature.shape)
            oldCursor = arcpy.SearchCursor(oriLayer)
            inAreaLEI=0
            error = 0
            for oldFeature in oldCursor:
                try:
                    insideArea=templeFeature.shape.intersect(oldFeature.shape,4)
                    inAreaLEI+=insideArea.area
                except Exception as e:
                    error = 1
                    print(e.message)
            if error == 0:
                leiDict[templeFeature.ID]=inAreaLEI/templeFeature.shape.area
            else:
                leiDict[templeFeature.ID]=999
        del tCursor
        del templeFeature
        del oldCursor
        del oldFeature

        arcpy.AddMessage('Step 4/4. Set LEI to new layer')
        print('Step 4/4. Set LEI to new layer')

        #Set LEI to the new layer
        newCursor = arcpy.UpdateCursor(newLayer)
        for newFeature in newCursor:
            if ".gdb" in iWorkspace:
                lei = leiDict[newFeature.OBJECTID]
            else:
                lei = leiDict[newFeature.FID]
            newFeature.LEI=lei*100
            newCursor.updateRow(newFeature)
        del newCursor
        del newFeature
    
        arcpy.AddMessage('Finished!')
        print('Finished!')

    except Exception as e:
        print(e.message)
        arcpy.AddMessage("Error: " + str(e.message))

		
		
if __name__=="__main__":

    """
    iWorkspace="F:\Projects\Liu\LEI\LEI\Data\Exercise1"
    oriShp="Old.shp"
    newShp="New.shp"
    dist=200
    """

    LEI(iWorkspace, oriShp, newShp, dist)
