# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Expandability
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import arcpy, sys

arcpy.env.overwriteOutput=1
built_area=sys.argv[1]
unbuildable_areas=sys.argv[2]
out_dir=sys.argv[3]
percision=float(sys.argv[4])

input_name=arcpy.Describe(built_area).name

if input_name.endswith('.shp'):
    input_name=input_name[:-4]

arcpy.AddMessage(out_dir+"/{0}_dissolved_buffer".format(input_name))
dsc=arcpy.Describe(out_dir)
if dsc.dataType=="Workspace":
    buffer_area=out_dir+"/buffer_Area"
    dissolved_buffer=out_dir+"/{0}_dissolved_buffer".format(input_name)
    unbuildable_merged=out_dir+"/{0}_unbuildable_merged".format(input_name)
    buildable_buffer=out_dir+"/{0}_buildable_buffer".format(input_name)
elif dsc.dataType=="Folder":
    buffer_area=out_dir+"/buffer_Area.shp"
    dissolved_buffer=out_dir+"/{0}_dissolved_buffer.shp".format(input_name)
    unbuildable_merged=out_dir+"/{0}_unbuildable_merged.shp".format(input_name)
    buildable_buffer=out_dir+"/{0}_buildable_buffer.shp".format(input_name)

sp=arcpy.Describe(built_area).spatialReference
linear_unit=sp.linearUnitName.lower()
fields=arcpy.ListFields(built_area)
fields=[i.name for i in fields]

if "AREA" not in fields:
    arcpy.AddField_management (built_area, "AREA","DOUBLE")
if "Perimeter" not in fields:
    arcpy.AddField_management (built_area, "Perimeter","DOUBLE")
if "OffsetDist" not in fields:
    arcpy.AddField_management (built_area, "OffsetDist","DOUBLE")
if "ExpanRate" not in fields:
    arcpy.AddField_management (built_area, "ExpanRate","DOUBLE")

arcpy.CalculateField_management (built_area, "AREA",'!shape.area@square{0}!'.format(linear_unit),'PYTHON')
arcpy.CalculateField_management (built_area, "Perimeter",'!shape.length@{0}!'.format(linear_unit),'PYTHON')
arcpy.CalculateField_management (built_area, "OffsetDist",'!AREA!/!Perimeter!','PYTHON')
                                 
arcpy.Buffer_analysis (built_area, buffer_area, "OffsetDist", "OUTSIDE_ONLY","FLAT","NONE")
arcpy.CalculateField_management (buffer_area, "ExpanRate",\
                                 '!shape.area@square{0}!/!AREA!'.format(linear_unit),'PYTHON')

criteria=0
while criteria==0:
    cursor_in_main=arcpy.da.UpdateCursor(built_area,['AREA','ExpanRate'])
    cursor_in_buffer=arcpy.da.SearchCursor(buffer_area,['ExpanRate'])
    for row_main in cursor_in_main:
        row_buffer=cursor_in_buffer.next()
        row_main[1]=row_buffer[0]
        cursor_in_main.updateRow(row_main)
    del cursor_in_main
    del cursor_in_buffer
    crit=1
    a=1
    cursor=arcpy.da.SearchCursor(built_area,['ExpanRate'])
    for row in cursor:
        if row[0]>(1+percision) or row[0]<(1-percision):
            a=0
        crit=crit*a
    del cursor
    criteria=crit
    arcpy.CalculateField_management (built_area, "OffsetDist",'!OffsetDist!/!ExpanRate!','PYTHON')
    arcpy.Buffer_analysis (built_area, buffer_area, "OffsetDist", "OUTSIDE_ONLY","FLAT","NONE")
    arcpy.CalculateField_management (buffer_area, "ExpanRate",\
                                     '!shape.area@square{0}!/!AREA!'.format(linear_unit),'PYTHON')
        

arcpy.Dissolve_management(buffer_area,dissolved_buffer,"#","#","MULTI_PART")

if unbuildable_areas!='false' and unbuildable_areas!='#':

    arcpy.Merge_management(unbuildable_areas, unbuildable_merged)
    arcpy.Erase_analysis(dissolved_buffer,unbuildable_merged,buildable_buffer)

    if "AREA" not in fields:
        arcpy.AddField_management (buildable_buffer, "AREA","DOUBLE")

    arcpy.CalculateField_management (buildable_buffer, "AREA",'!shape.area@square{0}!'.format(linear_unit),'PYTHON')

    cursor=arcpy.da.SearchCursor(buildable_buffer,['AREA'])
    row=cursor.next()
    buildable_A=row[0]
    del cursor

    built_A=0
    cursor=arcpy.da.SearchCursor(built_area,['AREA'])
    for row in cursor:
        built_A=built_A+row[0]
    del cursor

    expandability=buildable_A/built_A

    arcpy.AddMessage("Expandability: "+str(expandability))

    print expandability

    arcpy.SetParameterAsText(4,buildable_buffer)
    
else:
    arcpy.AddMessage("Expandability: 1")
    arcpy.SetParameterAsText(4,dissolved_buffer)
