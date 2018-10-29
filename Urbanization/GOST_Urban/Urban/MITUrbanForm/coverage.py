# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Coverage
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import sys, arcpy
arcpy.env.overwriteOutput=1
def getarea(fc):
    fields=arcpy.ListFields(fc)
    if "AREA" not in fields:
        arcpy.AddField_management (fc, "AREA","DOUBLE")
    arcpy.CalculateField_management (fc, "AREA", '!shape.area@squarekilometers!','PYTHON')
    fc_area=0
    cursor = arcpy.SearchCursor(fc)
    for row in cursor:
        fc_area=fc_area+row.getValue('AREA')
    del cursor
    return fc_area


blt=sys.argv[1]
unb=sys.argv[2]
out_dir=sys.argv[3]

dsc=arcpy.Describe(out_dir)

if dsc.dataType=="Workspace":
    convex=out_dir+"/"+"convex"
    unb_merged=out_dir+"/"+"unb_merge"
    unb_clip=out_dir+"/"+"unb_clip"
    
elif dsc.dataType=="Folder":
    convex=out_dir+"/"+"convex.shp"
    unb_merged=out_dir+"/"+"unb_merge.shp"
    unb_clip=out_dir+"/"+"unb_clip.shp"
    
arcpy.MinimumBoundingGeometry_management (blt, convex, 'CONVEX_HULL', 'ALL')
convex_area=getarea(convex)
blt_area=getarea(blt)

if unb!='false' and unb!='#':
    arcpy.AddMessage(str(unb))

    arcpy.Merge_management(unb, unb_merged)
    arcpy.Clip_analysis (unb_merged, convex,unb_clip)
    
    
    unb_area=getarea(unb_clip)
    msg=str(blt_area/(convex_area-unb_area))
    arcpy.AddMessage('Coverage: '+msg)
    arcpy.Delete_management(unb_merged)

    arcpy.SetParameterAsText(3,convex)
    arcpy.SetParameterAsText(4,unb_clip)
    
else:
    msg=str(blt_area/convex_area)
    arcpy.AddMessage('Coverage: '+msg)
    arcpy.SetParameterAsText(3,convex)
    
print msg