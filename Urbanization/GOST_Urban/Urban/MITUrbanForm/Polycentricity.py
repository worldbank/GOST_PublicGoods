# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Polycentricity
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import arcpy, sys,math, os, numpy
arcpy.env.overwriteOutput=1
in_data=sys.argv[1]
population_field=sys.argv[2]
unit=sys.argv[3]
out_dir=sys.argv[4]
dsc=arcpy.Describe(in_data)

dsc_directory=arcpy.Describe(out_dir)
if dsc_directory.dataType=="Workspace":
    centers=out_dir+"/centers"
    polygon_above_mean=out_dir+"/ply_abvMean"
    polygon_join_weight=out_dir+"/ply_join"
    polygon_full=out_dir+"/ply_abvMn_f"
    polygon_dissolved=out_dir+"/ply_dissolved"
    fc_convex=out_dir+"/convex"
    

    
elif dsc_directory.dataType=="Folder":
    centers=out_dir+"/centers.shp"
    polygon_above_mean=out_dir+"/ply_abvMean.shp"
    polygon_join_weight=out_dir+"/ply_join.shp"
    polygon_full=out_dir+"/ply_abvMn_f.shp"
    polygon_dissolved=out_dir+"/ply_dissolved.shp"
    fc_convex=out_dir+"/convex.shp"


total_weight=0
cursor = arcpy.SearchCursor(in_data)
for row in cursor:
    total_weight=total_weight+row.getValue(population_field)
del cursor

threshold_ratio=10/(total_weight)**0.5
threshold_weight=threshold_ratio*total_weight

print total_weight
if dsc.shapeType=="Point":
    search_radius=float(sys.argv[5])
    #raster cell size
    cell_size=sys.argv[6]
    extent=sys.argv[7]
    extracted_raster=out_dir+"/ex_raster"
    raster=out_dir+"/rastDen"

    if extent=='' or extent=='#':
        if dsc_directory.dataType=="Workspace":
            extent=out_dir+"/extent"
        elif dsc_directory.dataType=="Folder":
            extent=out_dir+"/extent.shp"
        arcpy.MinimumBoundingGeometry_management(in_data, fc_convex,'CONVEX_HULL')
        arcpy.Buffer_analysis (fc_convex, extent, search_radius)

    arcpy.CheckOutExtension("Spatial")
    
    #creat density raster (uses kernel density)
    outKDens=arcpy.sa.KernelDensity(in_data, population_field, cell_size, search_radius, unit)
    outKDens.save(raster)
    ext_rast_obj=arcpy.sa.ExtractByMask (raster, extent)
    ext_rast_obj.save(extracted_raster)
    #get the mean value (denisty) of raster
    
    meanObj=arcpy.GetRasterProperties_management(extracted_raster,"MEAN")
    mean=float(meanObj.getOutput(0))
    stdObj=arcpy.GetRasterProperties_management(extracted_raster,"STD")
    std=float(stdObj.getOutput(0))
    expression="'{0}'>{1}".format(extracted_raster,mean+2*std)

    #define the path-file-name for intermediate files
    raster_above_2stdFromMean=out_dir+"/rast_abv2std"

    #get the-above-mean parts of the density raster (creates a binary raster above=1 below=0)
    arcpy.gp.RasterCalculator(expression, raster_above_2stdFromMean)

    #define the path-file-name for intermediate files
    
    arcpy.RasterToPolygon_conversion (raster_above_2stdFromMean, polygon_full)
    
    print mean    
    arcpy.MakeFeatureLayer_management(polygon_full, "lyr", '"GRIDCODE"=1')
    arcpy.CopyFeatures_management("lyr", polygon_above_mean)
    field_map="""ID "ID" true true false 10 Double 0 10 ,First,#,\
    {0},ID,-1,-1;{1} "{1}" true true false 10 Double 0 10 ,\
    Sum,#,{2},{1},-1,-1""".format(polygon_above_mean,population_field,in_data)

    arcpy.SpatialJoin_analysis (polygon_above_mean, in_data, polygon_join_weight,"JOIN_ONE_TO_ONE",\
                          "KEEP_ALL",field_map,"INTERSECT")
    
    arcpy.Delete_management(polygon_full)
    arcpy.Delete_management(raster_above_2stdFromMean)

elif dsc.shapeType=="Polygon":
    fields=arcpy.ListFields(in_data)
    if "AREA_D" not in fields:
        arcpy.AddField_management (in_data, "AREA_D","DOUBLE")
    unit_fieldCalc=(unit.lower()).replace("_","")
    arcpy.CalculateField_management (in_data, "AREA_D",'!shape.area@{0}!'.format(unit_fieldCalc),'PYTHON')
    if "DENSITY" not in fields:
        arcpy.AddField_management (in_data, "DENSITY","DOUBLE")
    arcpy.CalculateField_management (in_data, "DENSITY", '!{0}!/!AREA_D!'.format(population_field),'PYTHON')
    sum_density=0
    count_dn=0
    cursor = arcpy.SearchCursor(in_data)
    for row in cursor:
        count_dn=count_dn+1
        sum_density=sum_density+row.getValue("DENSITY")
    del cursor
    mean_density=sum_density/count_dn
    
    cursor = arcpy.SearchCursor(in_data)
    # standard deviation
    sum_of_errors=0
    for row in cursor:
     sum_of_errors=sum_of_errors+((row.getValue("DENSITY")-mean_density)**2)
    del cursor
    std=(sum_of_errors/count_dn)**0.5
    arcpy.MakeFeatureLayer_management(in_data,"lyr",'"DENSITY">{0}'.format(mean_density+(2*std)))
    arcpy.CopyFeatures_management("lyr", polygon_above_mean)

    
    arcpy.Dissolve_management(polygon_above_mean,polygon_dissolved,"#","#",\
                              "SINGLE_PART","DISSOLVE_LINES")
    
    field_map="""ID "ID" true true false 10 Double 0 10 ,First,#,\
    {0},ID,-1,-1;{1} "{1}" true true false 10 Double 0 10 ,\
    Sum,#,{2},{1},-1,-1""".format(polygon_dissolved,population_field,in_data)
    
    arcpy.SpatialJoin_analysis (polygon_dissolved, in_data, polygon_join_weight,"JOIN_ONE_TO_ONE",\
                                "KEEP_ALL",field_map,"CONTAINS")


    arcpy.Delete_management(polygon_dissolved)
    



arcpy.MakeFeatureLayer_management(polygon_join_weight, "lyr_2",\
                                  '"{0}">{1}'.format(population_field,threshold_weight))
arcpy.CopyFeatures_management("lyr_2", centers)
    
entropy=0
result=arcpy.GetCount_management(centers)
count=int(result.getOutput(0))
##field = arcpy.da.TableToNumPyArray (centers, population_field, skip_nulls=True)
##centers_total_weight = field[population_field].sum()
centers_total_weight=0
cursor=arcpy.SearchCursor(centers)
for row in cursor:
    centers_total_weight=centers_total_weight+row.getValue(population_field)
del cursor

print count

if count==1:
    arcpy.AddMessage("case study is monocentric")
    arcpy.SetParameterAsText(6,centers)

elif count==0:
    arcpy.AddMessage("no center could be detected  by the default criteria")

else:
    center_jobs=0
    #computing entropy (zi*ln(zi))/ln(N)
    cursor = arcpy.SearchCursor(centers)
    for row in cursor:
        Z=row.getValue(population_field)/float(centers_total_weight)
        entropy=entropy+(((Z)*math.log(Z))/(math.log(count)))
    del cursor

    #computing entropy PC=(EI*N*Rc)
    print -entropy
    PC=-entropy*count*(centers_total_weight/float(total_weight))

    arcpy.AddMessage("Polycentricity= "+str(PC))
    arcpy.SetParameterAsText(7,centers)

arcpy.Delete_management(polygon_join_weight)
arcpy.Delete_management(polygon_above_mean)