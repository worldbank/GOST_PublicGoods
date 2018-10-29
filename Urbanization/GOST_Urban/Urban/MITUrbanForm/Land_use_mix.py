# ------------------------------------------------------------------------------
# Metropolitan Form Analysis Toolbox-Land-use Mix
# Credit: Reza Amindarbari, Andres Sevtsuk
# City Form Lab
# for more infomation see:
# Amindarbari, R., Sevtsuk, A., 2013, "Measuring Growth and Change in Metropolitan Form"
# presented at the Urban Affairs Association annual conference in San Francisco
# ------------------------------------------------------------------------------
import arcpy,sys,os
arcpy.env.overwriteOutput=1
def pointGrid(fc,outDir,unitCoversionRatio,gridSize):
    
    #dissolved_urban_extent=outDir+"/dissolved_urban_extent.shp"
    #arcpy.Dissolve_management(fc,dissolved_urban_extent,"#","#","MULTI_PART")
    fileName=os.path.basename(fc)
    dscObj=arcpy.Describe(fc)
    refExtent=dscObj.extent

    Y_Max=refExtent.YMax
    X_Max=refExtent.XMax
    Y_Min=refExtent.YMin
    X_Min=refExtent.XMin

    Y_length=(Y_Max-Y_Min)
    X_length=(X_Max-X_Min)

    No_Cols=X_length/gridSize
    No_Rows=Y_length/gridSize
    dsc=arcpy.Describe(outDir)
    if dsc.dataType=="Workspace":
        fishnet=outDir+"/fishnet"
        fishnet_label=outDir+"/fishnet_label"
        sample_point=outDir+"/{0}_sample_point".format(fileName)
    elif dsc.dataType=="Folder":
        fishnet=outDir+"/fishnet.shp"
        fishnet_label=outDir+"/fishnet_label.shp"
        sample_point=outDir+"/{0}_sample_point.shp".format(fileName)
        
    arcpy.CreateFishnet_management(fishnet,"{0} {1}".format(X_Min,Y_Min),\
                                   "{0} {1}".format(X_Min,Y_Max),(gridSize),(gridSize),\
                                   No_Rows,No_Cols,"#","LABELS","#","POLYGON")

    spReference=dscObj.spatialReference
    arcpy.DefineProjection_management(fishnet_label,spReference)
    arcpy.Clip_analysis (fishnet_label, fc, sample_point)
    

def MXcalculation(MTL,CW_LU_ratio_List,Fix_LU_ratio_List,\
                  Var_LU_ratio_Dic,LU_list,LU_area_local_List,Total_Area_Local):
    MXReport=''
    for i in MTL:
        
        if i==1:
            List=CW_LU_ratio_List
            
        if i==2:
            List=Fix_LU_ratio_List
            
        if i==3:
            indexMax=LU_area_local_List.index(max(LU_area_local_List))
            land_use_max=LU_list[indexMax]
            List=Var_LU_ratio_Dic[land_use_max]
            
        if i!=0:
            Ai=0 #sum_areas_of_interest
            counter=0
            Mi=1 #Matchin index
            if sum(LU_area_local_List)==0:
                MX_local=-9.9999
            else:
                for area in LU_area_local_List:
                    Ai=Ai+area
                    Mi=Mi*(1-abs(List[counter]-(area/sum(LU_area_local_List))))
                    counter=counter+1
                     
                Si=Ai/Total_Area_Local #Share of uses of interest in total
                MX_local=Si*Mi #local land-use mix index
                
            MXReport=MXReport+","+str(MX_local)
            #print [(Ai/Total_Area_Local),MX_local,Mi]
    return MXReport

##############################################################################
in_data=sys.argv[1]
land_use_field=sys.argv[2]
out_dir=sys.argv[3]
target_land_uses=sys.argv[4]
unit=sys.argv[5]
buffer_area=float(sys.argv[6]) #typically 1 sq-kilometer

if sys.argv[7]=='true':
    city_wide_land_use_distribution=sys.argv[7] #Boolean
else:
    city_wide_land_use_distribution=""
    
if sys.argv[8]=='true':
    fixed_expected_land_use_distribution=sys.argv[8] #Boolean
else:
    fixed_expected_land_use_distribution=""

if sys.argv[9]=='true':
    varying_expected_land_use_distribution=sys.argv[9] #Boolean
else:
    varying_expected_land_use_distribution=""
    
#MoranS_I=sys.argv[9] #Boolean


name=os.path.basename(in_data)
pi=3.14159265359

conversionDic={'Mile':{'Mile':1,'Meter':1609.34,'Foot':5280,'Kilometer':1.60934},\
               'Meter':{'Mile':0.000621371,'Meter':1,'Foot':3.28084,'Kilometer':0.001},\
               'Foot':{'Mile':0.000189394,'Meter':0.3048,'Foot':1,'Kilometer':0.0003048},\
               'Kilometer':{'Mile':0.621371,'Meter':1000,'Foot':3280.84,'Kilometer':1}}

if os.path.exists(out_dir+'/schema.ini'):
    os.remove(out_dir+'/schema.ini')
sp=arcpy.Describe(in_data).spatialReference
linear_unit=sp.linearUnitName
if linear_unit=='Foot_US':
    linear_unit='Foot'
    
conversionRatio= conversionDic[linear_unit][unit]

resolution_in_user_unit=(buffer_area/pi)**0.5 #distance between points (it is also the analysis radius)
resolution=resolution_in_user_unit/conversionRatio  #in map unit

target_land_uses_list=target_land_uses.split(';')

fields=arcpy.ListFields(in_data)
if "AREA" not in fields:
    arcpy.AddField_management (in_data, "AREA_LU","DOUBLE")
arcpy.CalculateField_management (in_data, "AREA_LU","!shape.area@square{0}!".format(unit.lower()),'PYTHON')

MX_report_header=''
measurement_type_list=[0,0,0]
sum_MX_locals=[]
MX_undefined=[]
messageTitles=[]

city_wide_land_use_ratio_List=[]
if city_wide_land_use_distribution:
    land_use_areas_List=[]
    messageTitles.append('Avergae land-use mix(ref. city-wide land-use ratios):')
    for i in range(1,len(target_land_uses_list)+1):
        land_use_areas_List.append(0)
    cursor=arcpy.da.SearchCursor(in_data,[land_use_field,"AREA_LU"])
    total_area=0
    for row in cursor:
        total_area=total_area+row[1]
        if row[0] in target_land_uses_list:
            land_use_areas_List[target_land_uses_list.index(row[0])]=land_use_areas_List[target_land_uses_list.index(row[0])]+row[1]

    city_wide_land_use_ratio_List=[(i/total_area) for i in land_use_areas_List]
    land_use_expected_CW=zip(target_land_uses_list,land_use_areas_List,city_wide_land_use_ratio_List)
    measurement_type_list[0]=1
    MX_report_header=MX_report_header+',MX_CW'
    sum_MX_locals.append(0)
    MX_undefined.append(0)
    
fixed_land_use_ratio_List=[]
if fixed_expected_land_use_distribution:
    fixed_land_use_ratio_str=sys.argv[10]
    fixed_land_use_ratio_List=fixed_land_use_ratio_str.split(';')
    fixed_land_use_ratio_List=[float(i) for i in fixed_land_use_ratio_List]
    
##    fixed_land_use_ratio=sys.argv[11]
##    land_use_areas_List_NA=[]
##    for i in range(1,len(target_land_uses_list)+1):
##        land_use_areas_List_NA.append('NA')
##    land_use_expected_FIXED=zip(target_land_uses_list,land_use_areas_List_NA,fixed_land_use_ratio_List)
    measurement_type_list[1]=2
    MX_report_header=MX_report_header+',MX_Fixed'
    sum_MX_locals.append(0)
    MX_undefined.append(0)
    messageTitles.append('Avergae land-use mix(ref. fixed ratios given by the user):')

varying_land_use_ratio_dic={}
if varying_expected_land_use_distribution:
    sum_MX_locals.append(0)
    MX_undefined.append(0)
    messageTitles.append('Avergae land-use mix(ref. varying ratios given by the user):')
    land_use_ratios_txt_file=sys.argv[11]
    file=open(land_use_ratios_txt_file,'r')
    varying_land_use_ratio_dic={}
    for line in file:
        if len(line.split('||'))>1:
            List=line.split('||')[1].split(';')
            #print List
            inner_list=[]
            main_key=line.split('||')[0]
            for i in List:
                ratio=i.split(':')[1].replace('\n','')
                inner_list.append(float(ratio))
                #print inner_dics
            varying_land_use_ratio_dic[main_key]=inner_list

##    land_use_areas_List_NA=[]
##    for i in range(1,len(target_land_uses_list)+1):
##        land_use_areas_List_NA.append('NA')
##    land_use_expected_FIXED=zip(target_land_uses_list,land_use_areas_List_NA,varying_land_use_ratio_dic)
    measurement_type_list[2]=3
    MX_report_header=MX_report_header+',MX_Varying'
    
pointGrid(in_data,out_dir,conversionRatio,resolution)

dsc=arcpy.Describe(out_dir)
if dsc.dataType=="Workspace":
    fishnet=out_dir+"/fishnet"
    fishnet_label=out_dir+"/fishnet_label"
    text_dir=os.path.dirname(out_dir)
    Land_Use_Dissolved=out_dir+"/Land_Use_Dissolved"
    sample_points=out_dir+"/{0}_sample_point".format(name)
    analysis_areas=out_dir+"/analysis_areas"
    intersects=out_dir+"/intersects"
    intersects_dissolved=out_dir+"/intersects_dissolved"
    summaryCSV=text_dir+'/summary.csv'
    schemaFile=text_dir+'/schema.ini'
    summaryDBF='{0}_summaryDBF'.format(name)
    ID="OBJECTID"
    intersect_FID="FID_analysis_areas"

if dsc.dataType=="Folder":
    fishnet=out_dir+"/fishnet.shp"
    fishnet_label=out_dir+"/fishnet_label.shp"
    Land_Use_Dissolved=out_dir+"/Land_Use_Dissolved.shp"
    sample_points=out_dir+"/{0}_sample_point.shp".format(name)
    analysis_areas=out_dir+"/analysis_areas.shp"
    intersects=out_dir+"/intersects.shp"
    intersects_dissolved=out_dir+"/intersects_dissolved.shp"
    summaryCSV=out_dir+'/summary.csv'
    schemaFile=out_dir+'/schema.ini'
    summaryDBF='{0}_summaryDBF.dbf'.format(name)
    ID="FID"
    intersect_FID="FID_analys"
    
arcpy.Dissolve_management(in_data,Land_Use_Dissolved,land_use_field,"#","MULTI_PART")
arcpy.Buffer_analysis(sample_points,analysis_areas,resolution)
#arcpy.Dissolve_management(in_data,Land_Use_Dissolved,land_use_field,"#","MULTI_PART")
arcpy.Intersect_analysis([analysis_areas,Land_Use_Dissolved],intersects)
arcpy.Dissolve_management(intersects,intersects_dissolved,[land_use_field,intersect_FID],"#","MULTI_PART")

fields=arcpy.ListFields(intersects_dissolved)
fields=[i.name for i in fields]
if "AREA_LU" not in fields:
    arcpy.AddField_management (intersects_dissolved, "AREA_LU","DOUBLE")
arcpy.CalculateField_management (intersects_dissolved, "AREA_LU",\
                                 "!shape.area@square{0}!".format(unit),'PYTHON')

infile = open(summaryCSV,'w')
target_land_uses_header=target_land_uses.replace(';',',')


infile.write("{0},Total_Area,".format(intersect_FID)+target_land_uses_header+MX_report_header+"\n")
cursor=sorted(arcpy.da.SearchCursor(intersects_dissolved,\
                                    [intersect_FID,land_use_field,'AREA_LU']))

a=0 
b=0
total_area=0
backup_area=0
land_use_areas_local_List=[]
for i in range(1,len(target_land_uses_list)+1):
    land_use_areas_local_List.append(0)
    
for row in cursor:
    a=a+row[0]
    b=b+1
    total_area=total_area+row[2]
    if (a/b)==row[0]:
        if row[1] in target_land_uses_list:
            #backup=land_use_areas_local_List[target_land_uses_list.index(row[1])]
            land_use_areas_local_List[target_land_uses_list.index(row[1])]=row[2]
            #print land_use_areas_local_List
    else:
        a=row[0]
        b=1
        total_area=total_area-row[2]
        #land_use_areas_local_List[target_land_uses_list.index(row[1])]=backup
        land_use_areas_local_List_str=[str(i) for i in land_use_areas_local_List]
        land_use_areas_local_str=','.join(land_use_areas_local_List_str)

        MX_report=MXcalculation(measurement_type_list,city_wide_land_use_ratio_List,fixed_land_use_ratio_List,\
              varying_land_use_ratio_dic,target_land_uses_list,\
                                land_use_areas_local_List,total_area)
        MX_report_list=MX_report.split(',')
        MX_report_list.pop(0)
        for i in range(len(MX_report_list)):
            if MX_report_list[i]=='-9.9999':
                MX_undefined[i]=MX_undefined[i]-9.9999

        sum_MX_locals=[x+float(y) for x,y in zip(sum_MX_locals,MX_report_list)]
        
    
        infile.write(str(row[0]-1)+","+str(total_area)+","+land_use_areas_local_str+MX_report+"\n")
        total_area=row[2]
        backup_area=row[2]
        backup_land_use=row[1]
        land_use_areas_local_List=[]
        for i in range(1,len(target_land_uses_list)+1):
            land_use_areas_local_List.append(0)
        if backup_area and (backup_land_use in target_land_uses_list):
            land_use_areas_local_List[target_land_uses_list.index(backup_land_use)]=backup_area
            backup_area=0

MX_report=MXcalculation(measurement_type_list,city_wide_land_use_ratio_List,fixed_land_use_ratio_List,\
                        varying_land_use_ratio_dic,target_land_uses_list,\
                        land_use_areas_local_List,total_area)

land_use_areas_local_List_str=[str(i) for i in land_use_areas_local_List]
land_use_areas_local_str=','.join(land_use_areas_local_List_str)
infile.write(str(row[0])+","+str(total_area)+","+land_use_areas_local_str+MX_report+"\n")
del cursor
del infile

no_Fields=len(target_land_uses_list)+3

##if os.path.exists(schemaFile):
##    os.remove(schemaFile)

schema=open(schemaFile,'w')
schema.write('[{0}]\n'.format(summaryCSV))
if city_wide_land_use_distribution:
    schema.write('Col{0}=MX_CW Double\n'.format(no_Fields))
if fixed_expected_land_use_distribution:
    if city_wide_land_use_distribution:
        schema.write('Col{0}=MX_Fixed Double\n'.format(no_Fields+1))
    else:
        schema.write('Col{0}=MX_Fixed Double\n'.format(no_Fields))
if varying_expected_land_use_distribution:
    if city_wide_land_use_distribution and fixed_expected_land_use_distribution:
        schema.write('Col{0}=MX_Varying Double\n'.format(no_Fields+2))
    elif (city_wide_land_use_distribution or fixed_expected_land_use_distribution)\
         and not(city_wide_land_use_distribution and fixed_expected_land_use_distribution):
        schema.write('Col{0}=MX_Varying Double\n'.format(no_Fields+1))
    else:
        schema.write('Col{0}=MX_Varying Double\n'.format(no_Fields))
del schema

arcpy.TableToTable_conversion(summaryCSV,out_dir,summaryDBF)
arcpy.JoinField_management (sample_points, ID,out_dir+"/"+summaryDBF,intersect_FID)

MX_undefined=[-i for i in MX_undefined]
sum_MX_locals=[x+y for x,y in zip(MX_undefined,sum_MX_locals)]
result=arcpy.GetCount_management(sample_points)
count=int(result.getOutput(0))

ave_MX=[x/(count-(y/9)) for x,y in zip(sum_MX_locals,MX_undefined)]
for title,average_MX in zip(messageTitles,ave_MX):
    arcpy.AddMessage(title)    
    arcpy.AddMessage(str(average_MX))

params = arcpy.GetParameterInfo()

scriptDir = os.path.dirname(sys.argv[0])
MX_CW = os.path.join(scriptDir, "Symbology_MX_CW.lyr")
MX_Fixed = os.path.join(scriptDir, "Symbology_MX_Fixed.lyr")
MX_Varying = os.path.join(scriptDir,"Symbology_MX_Varying.lyr")

if varying_expected_land_use_distribution:
    params[11].symbology=MX_Varying
elif fixed_expected_land_use_distribution:
    params[11].symbology=MX_Fixed
else:
    params[11].symbology=MX_CW


arcpy.SetParameterAsText(11,sample_points)
arcpy.Delete_management(fishnet)
arcpy.Delete_management(fishnet_label)
arcpy.Delete_management(Land_Use_Dissolved)
arcpy.Delete_management(analysis_areas)
arcpy.Delete_management(intersects)