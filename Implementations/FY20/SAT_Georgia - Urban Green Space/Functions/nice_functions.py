from functools import partial
import pyproj

from shapely.geometry import box, shape
from shapely.ops import transform

import pandas as pd

import gbdxtools
from gbdxtools import Interface
gbdx = Interface()

from collections import namedtuple


def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')

def largest(features, n=100):
    areas = []
    for f in features:
        s = shape(f['geometry'])
        areas.append((s.area, f))

    topN = [r for r in reversed(sorted(areas, key=lambda k: k[0]))][0:n]
    return zip(*topN)


def find(bbox, query):
    wkt = box(*bbox).wkt
    return gbdx.vectors.query(wkt, query, index=None)

def create_mask_from_vector(vector_data_path, cols, rows, geo_transform,
                            projection, target_value=1):
    """Rasterize the given vector (wrapper for gdal.RasterizeLayer)."""
    data_source = gdal.OpenEx(vector_data_path, gdal.OF_VECTOR)
    layer = data_source.GetLayer(0)
    driver = gdal.GetDriverByName('MEM')  # In memory dataset
    target_ds = driver.Create('', cols, rows, 1, gdal.GDT_UInt16)
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(projection)
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[target_value])
    return target_ds


def vectors_to_raster(file_paths, rows, cols, geo_transform, projection):
    """Rasterize the vectors in the given directory in a single image."""
    labeled_pixels = np.zeros((rows, cols))
    for i, path in enumerate(file_paths):
        label = i+1
        ds = create_mask_from_vector(path, cols, rows, geo_transform,
                                     projection, target_value=label)
        band = ds.GetRasterBand(1)
        labeled_pixels += band.ReadAsArray()
        ds = None
    return labeled_pixelsh


def get_projections_and_UTM(bbox_city):

    # convert area to string for use in request
    bbox_large_area_str = listToStringWithoutBrackets(bbox_city)

    # UTM zone and EPSG code calculator
    zone_cal = round((183+bbox_city[0])/6,0)
    EPSG = 32700-round((45+bbox_city[1])/90,0)*100+round((183+bbox_city[0])/6,0)

    UTM_EPSG_code ='EPSG:%i'  % (EPSG)

    # create reprojection variables
    project_utm = partial(
        pyproj.transform,
        pyproj.Proj(init='epsg:4326'), # source coordinate system
        pyproj.Proj(init=UTM_EPSG_code)) # destination coordinate system

    # create reprojection variables
    project_wgs = partial(
        pyproj.transform,
        pyproj.Proj(init=UTM_EPSG_code), # destination coordinate system
        pyproj.Proj(init='epsg:4326')) # source coordinate system


    return UTM_EPSG_code,project_utm,project_wgs


bbox_amsterdam = [4.7288558, 52.27813889999999, 5.0683903, 52.4311573]


def get_OSM_polygons(bbox_city,type_query = 'park',count = 1e4):

    if bbox_city == 0:
        return 'no bbox'



    # convert the bounding box into a well-known text format
    bbox_wkt = box(*bbox_city).wkt

    # for water use: AND attributes.natural:water
    # for grass/forest use: AND attributes.landuse:grass or forest
    # for footway use: AND attributes.highway:footway
    # for park use: item_type:Park
    #item_type:Park
    if type_query == 'park':

        query = "ingest_source:OSM AND attributes.leisure:park AND geom_type:Polygon "

    elif type_query == 'grass':

        query = "ingest_source:OSM AND (attributes.landuse:grass OR attributes.natural:grassland) AND geom_type:Polygon "

    elif type_query == 'forest':

        query = "ingest_source:OSM AND (attributes.landuse:forest OR attributes.natural:wood) AND geom_type:Polygon "

    elif type_query == 'water':

        query = "ingest_source:OSM AND attributes.natural:water AND geom_type:Polygon "

    elif type_query == 'footway':

        query = "ingest_source:OSM AND attributes.highway:footway AND geom_type:Polygon "

    elif type_query == 'building':

        query = "ingest_source:OSM AND item_type:Building AND geom_type:Polygon "

    else :

        return 'wrong type_query'

    # query the gbdx.vectors

    geojson_obj = gbdx.vectors.query(bbox_wkt, query= query , index="vector-osm-*", count = count)

    return geojson_obj



def geojson_variable_check(geojson_obj,project_utm,min_size = 1):

    # create dataframe to save all data
    df = pd.DataFrame(columns=['id','OSM_id','item_type','name','geom_type','area','check'])

    # Convert geojson from OSM to shapely polygons
    geom_list = []
    for geojson in geojson_obj:
        geom = shape(geojson['geometry'])
        geom_list.append(geom)

    # Loop over all polygons and load information to dataframe
    i = 0;
    for r in geojson_obj:
        # get properties
        props = r['properties']
        # get geometry type
        geom_type = r['geometry']['type']


        # measure area if geometry is polygon otherwise set area to 0
        if geom_type == u'Polygon':

            utm_projected = transform(project_utm, geom_list[i])  # apply UTM projection

            area_obj = utm_projected.area/10000 # calculate area in ha
        else:
            area_obj = 0;

        # load all metadata to dataframe
        df = df.append({'id':i ,'OSM_id':props[u'id'],'item_type':props[u'item_type'][0],'name':props[u'name'],
                       'geom_type':geom_type,'area':area_obj,
                        'check':geom_type == u'Polygon' and area_obj > min_size},ignore_index=True) #set test for desirable parks

        i = i + 1




    # get indices of all parks that pass the test set above
    list_obj = df.loc[df.check == True]['id'];

    final_list = df.loc[df.check == True].reset_index() #.sort_values(by=['area'], ascending = False).reset_index()



    return final_list, geom_list

def image_query_check(bbox_wkt, park_utm, buffer_size, project_wgs,x_wgs,y_wgs):


    if (park_utm.buffer(buffer_size).exterior == None):

        message = "error"

        return message

    else:



        bbox_park_area_float = min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)

        # query images that cover the park
        # query = "(item_type:WV03_VNIR OR WV02 OR QB02 OR GE01)"

        query = "(item_type:WV03_VNIR OR WV02)"
        query += "AND NOT item_type:IDAHOImage AND item_type:DigitalGlobeProduct"
        results = gbdx.vectors.query(bbox_wkt, query)

        # calculate the overlap between image and park

        Rectangle = namedtuple('Rectangle', 'xmin ymin xmax ymax')

        ra = Rectangle(min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)) # set wgs park area in rectangle

        # function for calculating overlap
        def area_overlap(a, b):  # returns None if rectangles don't intersect
            dx = min(a.xmax, b.xmax) - max(a.xmin, b.xmin)
            dy = min(a.ymax, b.ymax) - max(a.ymin, b.ymin)
            if (dx>=0) and (dy>=0):
                return dx*dy




        # create dataframe for image metadata
        images_df = pd.DataFrame(columns=['id','month','year','type','cloud cover','overlap','check'])

        # loop over all found images and load metadata in dataframe
        # also perform a test to find the right images
        for r in results:
            props = r['properties']
            geom = r['geometry']

            coordinates_image = geom['coordinates'][0][0]
            df_overlap = pd.DataFrame.from_records(coordinates_image,columns = ['x','y'])

            rb = Rectangle(df_overlap.x.min(), df_overlap.y.min(), df_overlap.x.max(), df_overlap.y.max())

            total_area = ((ra.xmax-ra.xmin)*(ra.ymax-ra.ymin))

            total_overlap = area_overlap(ra, rb)

            if total_overlap is None:
                fraction_overlap = 0
            else:
                fraction_overlap = area_overlap(ra, rb)/total_area


            images_df = images_df.append({'id': props['attributes']['catalogID'],'month':props['item_date'][5:7],'year':props['item_date'][0:4],
                           'type':props['item_type'][1],
                           'cloud cover':props['attributes']['cloudCover_int'],'overlap': fraction_overlap,
                            # perform check on, cloud cover, month (summer), and overlap
                            'check':props['attributes']['cloudCover_int'] < 5 and int(props['item_date'][5:7]) >= 5 and int(props['item_date'][5:7]) <= 11 and
                           fraction_overlap >.9},ignore_index=True)


        # select only the images that survive the test
        selection_images = images_df.loc[images_df.check == True].reset_index()

        # group by image id because
        selection_images =  selection_images.groupby(['id'],sort=['year','month'],as_index=False).first().sort_values(['year','month'], ascending=False).reset_index()



        return selection_images
