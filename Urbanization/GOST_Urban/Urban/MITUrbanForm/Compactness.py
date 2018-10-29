# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Compactness
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import arcpy, sys, os, math
sys.path.append(r"C:\Users\wb411133\Box Sync\AAA_BPS\Code\GOST")
from GOSTRocks.misc import tPrint
from GOSTRocks.arcpyMisc import tryAddField

arcpy.env.overwriteOutput=1

def rawCompactness(fc, field, unitCoversionRation, beta):
    fields=arcpy.ListFields(fc)
    if "XCord" not in fields:
        arcpy.AddField_management (fc, "XCord","DOUBLE")
    if "YCord" not in fields:
        arcpy.AddField_management (fc, "YCord","DOUBLE")
    arcpy.CalculateField_management (fc, "XCord","float(!SHAPE.CENTROID!.split()[0])",'PYTHON')
    arcpy.CalculateField_management (fc, "YCord","float(!SHAPE.CENTROID!.split()[1])",'PYTHON')


    cursorO=arcpy.da.SearchCursor(fc,['XCord','YCord',field])
    cursorD=arcpy.da.SearchCursor(fc,['XCord','YCord',field])

    Sum_G_weighted=0
    Sum_population=0
    if beta == 0:
        for rowO in cursorO:
            G=0
            for rowD in cursorD:
                dist=unitCoversionRation*((((rowD[0]-rowO[0])**2)+((rowD[1]-rowO[1])**2))**0.5)
                if dist:
                    G=G+(int(rowD[2])/dist)
            Sum_G_weighted=Sum_G_weighted+(G*int(rowO[2]))
            Sum_population=Sum_population+(int(rowO[2]))
    else:
        for rowO in cursorO:
            G=0
            for rowD in cursorD:
                dist=unitCoversionRation*((((rowD[0]-rowO[0])**2)+((rowD[1]-rowO[1])**2))**0.5)
                if dist:
                    G=G+(int(rowD[2])/(Math.e**(float(beta)*dist)))
            Sum_G_weighted=Sum_G_weighted+(G*int(rowO[2]))
            Sum_population=Sum_population+int(rowO[2])
            
    del cursorO
    del cursorD

    compactness=Sum_G_weighted/Sum_population
    return [compactness,Sum_population]

def pointGrid(polygon,outDir,unitCoversionRation,gridSize,total_population,type):
    dscObj=arcpy.Describe(polygon)
    refExtent=dscObj.extent

    Y_Max=refExtent.YMax
    X_Max=refExtent.XMax
    Y_Min=refExtent.YMin
    X_Min=refExtent.XMin

    Y_length=unitCoversionRation*(Y_Max-Y_Min)
    X_length=unitCoversionRation*(X_Max-X_Min)

    No_Cols=X_length/gridSize
    No_Rows=Y_length/gridSize
    dsc=arcpy.Describe(outDir)
    if dsc.dataType=="Workspace":
        fishnet=outDir+"/fishnet"
        fishnet_label=outDir+"/fishnet_label"
        reference_point=outDir+"/reference_point_{0}".format(type)
    elif dsc.dataType=="Folder":
        fishnet=outDir+"/fishnet.shp"
        fishnet_label=outDir+"/fishnet_label.shp"
        reference_point=outDir+"/reference_point_{0}.shp".format(type)
        
    arcpy.CreateFishnet_management(fishnet,"{0} {1}".format(X_Min,Y_Min),\
                                   "{0} {1}".format(X_Min,Y_Max),(gridSize/unitCoversionRation),(gridSize/unitCoversionRation),\
                                   No_Rows,No_Cols,"#","LABELS","#","POLYGON")

    spReference=dscObj.spatialReference
    
    arcpy.DefineProjection_management(fishnet_label,spReference)
    arcpy.Clip_analysis (fishnet_label, polygon, reference_point)
    fields=arcpy.ListFields(reference_point)
    if "population" not in fields:
        arcpy.AddField_management (reference_point, "population","DOUBLE")

    number_of_points = int(arcpy.GetCount_management(reference_point).getOutput(0))
    arcpy.CalculateField_management (reference_point, "population",total_population/number_of_points,'PYTHON')
##############################################################################

def processUrbanExtents(extents, popRaster, tempFolder, fieldPref="c"):
    allCompactness = []
    filesToDelete = []
    #add two output fields
    tryAddField(extents, "%sCom" % fieldPref, "FLOAT")
    tryAddField(extents, "%sNCom" % fieldPref, "FLOAT")
    with arcpy.da.UpdateCursor(extents, ["OID@", "%sCom" % fieldPref, "%sNCom" % fieldPref]) as cursor:
        for featRow in cursor:
            tRaster = os.path.join(tempFolder, "pop_%s" % featRow[0])
            tRasterPoints = os.path.join(tempFolder, "pop_%s_pts.shp" % featRow[0])
            tLayer = "pop_%s" % featRow[0]
            filesToDelete.append(tRaster)
            filesToDelete.append(tRasterPoints)
            filesToDelete.append(tRasterPoints)
            
            arcpy.MakeFeatureLayer_management(extents, tLayer, '"FID" = %s' % featRow[0])
            #Use the current feature to extract the current raster
            arcpy.Clip_management(popRaster, '#', tRaster, tLayer, '0', 'ClippingGeometry', 'MAINTAIN_EXTENT')
            try:
                arcpy.RasterToPoint_conversion(tRaster, tRasterPoints)
                compactness = main(tRasterPoints, "GRID_CODE", tempFolder)
                featRow[1] = compactness[0]
                featRow[2] = compactness[1]
            except:
                tPrint("Something went wrong with feature %s" % featRow[0])
                featRow[1] = 0
                featRow[2] = 0
            cursor.updateRow(featRow)
    for f in filesToDelete:
        arcpy.Delete_management(f)
    return(allCompactness)

def main(in_data, population_field, out_dir,
        unit="Kilometer", beta=0, norm_by_reference="None", #None without_Geographic_Constraints with_Geographic_Constraints Both
        reference_density=300, unbuildable='', #Geographic_Constraints determines the maximum difference between the area of reference after removing unbuildable area and the area of perfect circular reference
        percision=0.5):
    e=2.718281828459045
    pi=3.14159265359
    conversionDic={'Mile':{'Mile':1,'Meter':1609.34,'Foot':5280,'Kilometer':1.60934},\
                   'Meter':{'Mile':0.000621371,'Meter':1,'Foot':3.28084,'Kilometer':0.001},\
                   'Foot':{'Mile':0.000189394,'Meter':0.3048,'Foot':1,'Kilometer':0.0003048},\
                   'Kilometer':{'Mile':0.621371,'Meter':1000,'Foot':3280.84,'Kilometer':1}}

    sp = arcpy.Describe(in_data).spatialReference
    linear_unit = sp.linearUnitName
    #print linear_unit
    conversionRatio = conversionDic[linear_unit][unit]
    #print conversionRatio
    if linear_unit=='Foot_US':
        linear_unit='Foot'
    case_rawCompactness = rawCompactness(in_data, population_field, conversionRatio, beta)
    total_population = case_rawCompactness[1]
    raw_Compactness = case_rawCompactness[0]
    normalizedCompactness = raw_Compactness/total_population
    if norm_by_reference == "None":
        return([raw_Compactness, normalizedCompactness])
    '''
    print "raw compactness:", case_rawCompactness[0]
    arcpy.AddMessage("raw compactness: "+str(case_rawCompactness[0]))
    print "compactness normalized by population:", case_rawCompactness[0]/total_population
    arcpy.AddMessage("compactness normalized by population:"+str(case_rawCompactness[0]/total_population))
    '''
    dsc=arcpy.Describe(out_dir)
    if dsc.dataType=="Workspace":
        fishnet=out_dir+"/fishnet"
        fishnet_label=out_dir+"/fishnet_label"
        center=out_dir+"/center"
        reference=out_dir+"/reference"
        reference_point_c=out_dir+"/reference_point_C"
        unbuildable_merged=out_dir+"/unbuildable_merged"
        reference_w_gConstraints_step1=out_dir+"/reference_w_gConstraints_step1"
        reference_w_gConstraints_step2=out_dir+"/reference_w_gConstraints_step2"
        reference_w_gConstraints=out_dir+"/reference_w_gConstraints"
        reference_point_CG=out_dir+"/reference_point_CG"

    elif dsc.dataType=="Folder":
        fishnet=out_dir+"/fishnet.shp"
        fishnet_label=out_dir+"/fishnet_label.shp"
        center=out_dir+"/center.shp"
        reference=out_dir+"/reference.shp"
        reference_point_c=out_dir+"/reference_point_C.shp"
        unbuildable_merged=out_dir+"/unbuildable_merged.shp"
        reference_w_gConstraints_step1=out_dir+"/reference_w_gConstraints_step1.shp"
        reference_w_gConstraints_step2=out_dir+"/reference_w_gConstraints_step2.shp"
        reference_w_gConstraints=out_dir+"/reference_w_gConstraints.shp"
        reference_point_CG=out_dir+"/reference_point_CG.shp"

    if norm_by_reference!="None":
        resolution=float(sys.argv[8])# for converting the reference polygon to point 
        referenceArea=total_population/float(reference_density)
        referenceRadius=((referenceArea/math.pi)**0.5)/conversionRatio
        arcpy.MeanCenter_stats(in_data,center, population_field)
        arcpy.Buffer_analysis(center,reference,referenceRadius)
        
        if norm_by_reference=="without_Geographic_Constraints" or norm_by_reference=="Both":
            pointGrid(reference,out_dir,conversionRatio,resolution,total_population,"C")
            reference_rawCompactness=rawCompactness(reference_point_c,"population",conversionRatio)
            print "raw compactness of the reference (without geographic constraints): ",reference_rawCompactness[0]
            arcpy.AddMessage("raw compactness of the reference (without geographic constraints): "+str(reference_rawCompactness[0]))
            #print reference_rawCompactness[1]
            print "compactness normalized by a circular reference:",(case_rawCompactness[0]/reference_rawCompactness[0])
            arcpy.AddMessage("compactness normalized by the reference(without geographic constraints): "+str(case_rawCompactness[0]/reference_rawCompactness[0]))
            arcpy.SetParameterAsText(10,reference)
            
        if norm_by_reference=="with_Geographic_Constraints" or norm_by_reference=="Both":
            arcpy.Merge_management (unbuildable, unbuildable_merged)
            arcpy.Erase_analysis (reference, unbuildable_merged, reference_w_gConstraints_step1)
            fields=arcpy.ListFields(reference_w_gConstraints_step1)
            if "AREA" not in fields:
                arcpy.AddField_management (reference_w_gConstraints_step1, "AREA","DOUBLE")
            arcpy.CalculateField_management (reference_w_gConstraints_step1, "AREA","!shape.area@square{0}!".format(unit),'PYTHON')
            cursor=arcpy.da.SearchCursor(reference_w_gConstraints_step1,"AREA")
            for row in cursor:
                referenceArea_w_gConstraints=row[0]
            del cursor
            
            buffer_dist=referenceRadius*(((referenceArea/referenceArea_w_gConstraints)**0.5)-1)
            arcpy.Buffer_analysis(reference_w_gConstraints_step1,reference_w_gConstraints_step2,buffer_dist) 

            expansionRate=0
            while expansionRate<(1-percision) or expansionRate>(1+percision):
                arcpy.Erase_analysis (reference_w_gConstraints_step2, unbuildable_merged, reference_w_gConstraints)
                fields=arcpy.ListFields(reference_w_gConstraints)
                if "AREA" not in fields:
                    arcpy.AddField_management (reference_w_gConstraints, "AREA","DOUBLE")            
                arcpy.CalculateField_management (reference_w_gConstraints, "AREA","!shape.area@square{0}!".format(unit),'PYTHON')
                cursor=arcpy.da.SearchCursor(reference_w_gConstraints,"AREA")
                for row in cursor:
                    referenceArea_w_gConstraints=row[0]
                del cursor
                expansionRate=referenceArea/referenceArea_w_gConstraints
                buffer_dist=buffer_dist*expansionRate
                arcpy.Buffer_analysis(reference_w_gConstraints_step1,reference_w_gConstraints_step2,buffer_dist)
                
            pointGrid(reference_w_gConstraints,out_dir,conversionRatio,resolution,total_population,"CG")
            reference_rawCompactness_w_gConstraints=rawCompactness(reference_point_CG,"population",conversionRatio)
            print "raw compactness of the reference (with geographic constraints)",reference_rawCompactness_w_gConstraints[0]
            arcpy.AddMessage("raw compactness of the reference (with geographic constraints): "+str(reference_rawCompactness_w_gConstraints[0]))
            #print reference_rawCompactness_w_gConstraints[1]
            print "compactness normalized by the reference(with geographic constraints):",(case_rawCompactness[0]/reference_rawCompactness_w_gConstraints[0])
            arcpy.AddMessage("compactness normalized by the reference(with geographic constraints):"+str(case_rawCompactness[0]/reference_rawCompactness_w_gConstraints[0]))
            arcpy.SetParameterAsText(11,reference_w_gConstraints)
            arcpy.Delete_management(fishnet_label)
            arcpy.Delete_management(fishnet)
            arcpy.Delete_management(center)
            arcpy.Delete_management(reference_w_gConstraints_step1)
            arcpy.Delete_management(reference_w_gConstraints_step2)
            arcpy.Delete_management(unbuildable_merged)
            