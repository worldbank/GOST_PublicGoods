### other libraries

import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
import math
import copy
from shapely.ops import transform
from shapely.geometry import mapping, Polygon, box, shape
from skimage.transform import resize
import matplotlib.pyplot as plt


from collections import namedtuple

from gbdxtools import Interface
from gbdxtools.task import env
from gbdxtools import CatalogImage

gbdx = Interface()

# convert a list of numbers to a string list without brackets or parentheses
def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')

#  This funtion resizes images with Skimage so they are the same size, i resamples the largest to a courser resolution. Also returns a number of the image that was reshaped because it will be transformed into an array and does not hold all information anymore.
def reziseimages(image1,image2):

    #first lets see which of the images is the largest (highest Resolution) and then resize that one
    if image1.shape[2] > image2.shape[2]:
        image1 = resize(image1.astype(np.float), image2.shape, order=1, mode='constant', clip=True)
        reshaped = 1
    else:
        image2 = resize(image2.astype(np.float), image1.shape, order=1, mode='constant', clip=True)
        reshaped = 2
    return image1,image2, reshaped

# input is a shapely polygon of the park of interest and output is two dataframes with listed images and the most recent image for Summer and Winter with least amount of cloud cover
def get_SumWin(park_shape,UTM_EPSG_code):

#     Getting the bounding boxes right
    x_wgs,y_wgs = park_shape.exterior.xy

    bbox_park_area = str([min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)])
    bbox_park_area_str = listToStringWithoutBrackets(bbox_park_area)

    bbox_park_area = min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)
    bbox_wkt = box(*bbox_park_area).wkt

    aoi = bbox_wkt

# Querying for WV03 or WV02 within the aoi
    query = "(item_type:WV03_VNIR OR WV02)"
    query += " AND NOT item_type:IDAHOImage AND item_type:DigitalGlobeProduct"
    results = gbdx.vectors.query(aoi, query)


    from collections import namedtuple
    Rectangle = namedtuple('Rectangle', 'xmin ymin xmax ymax')

    ra = Rectangle(min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs))

    # intersection here is (3, 3, 4, 3.5), or an area of 1*.5=.5

    def area(a, b):  # returns None if rectangles don't intersect
        dx = min(a.xmax, b.xmax) - max(a.xmin, b.xmin)
        dy = min(a.ymax, b.ymax) - max(a.ymin, b.ymin)
        if (dx>=0) and (dy>=0):
            return dx*dy

# set cloud check

    images_df = pd.DataFrame(columns=['id','month','year','type','resolution','cloud cover','overlap','check'])
    images_df2 = pd.DataFrame(columns=['id','month','year','type','resolution','cloud cover','overlap','check'])

    for r in results:
        props = r['properties']
        geom = r['geometry']

        coordinates_image = geom['coordinates'][0][0]
        df_overlap = pd.DataFrame.from_records(coordinates_image,columns = ['x','y'])

        rb = Rectangle(df_overlap.x.min(), df_overlap.y.min(), df_overlap.x.max(), df_overlap.y.max())

        total_area = ((ra.xmax-ra.xmin)*(ra.ymax-ra.ymin))

        total_overlap = area(ra, rb)

        if total_overlap is None:
            fraction_overlap = 0
        else:
            fraction_overlap = area(ra, rb)/total_area


        images_df = images_df.append({'id': props['attributes']['catalogID'],'month':props['item_date'][5:7],'year':props['item_date'][0:4],'type':props['item_type'][1],'resolution':props['attributes']['resolution_dbl'],'cloud cover':props['attributes']['cloudCover_int'],'overlap': fraction_overlap,'check':props['attributes']['cloudCover_int'] < 100 and int(props['item_date'][5:7]) >= 5 and int(props['item_date'][5:7]) <= 8 and fraction_overlap >.9},ignore_index=True)

        images_df2 = images_df2.append({'id': props['attributes']['catalogID'],'month':props['item_date'][5:7],'year':props['item_date'][0:4],'type':props['item_type'][1],'resolution':props['attributes']['resolution_dbl'],'cloud cover':props['attributes']['cloudCover_int'],'overlap': fraction_overlap,'check':props['attributes']['cloudCover_int'] < 100 and int(props['item_date'][5:7]) <= 2 or int(props['item_date'][5:7]) >= 11 and fraction_overlap >.9},ignore_index=True)


    selections = images_df.loc[images_df.check == True].reset_index()
    selectionw = images_df2.loc[images_df2.check == True].reset_index()

    selections = selections.groupby(['id'],sort=['year','month','cloud cover'],as_index=False).first().sort_values(['year','month'], ascending=False).reset_index()
    del selections['level_0'], selections['index']

    selectionw = selectionw.groupby(['id'],sort=['year','month','cloud cover'],as_index=False).first().sort_values(['year','month'], ascending=False).reset_index()
    del selectionw['level_0'], selectionw['index']

# First check if images are not defective and then fetch them

    bbox = env.inputs.get('bbox', bbox_park_area_str)
#     Select the first image
    imagestatus ='unknown'
    iS=0
    iW=0



    if (selections.empty | selectionw.empty):


        error = 0

        return error, error, error, error


    else:




        while imagestatus =='unknown':
            # set catalog id from selection
            catalog_id_s = env.inputs.get('catalog_id', selections.id[iS])
            catalog_id_w = env.inputs.get('catalog_id', selectionw.id[iW])

        # collect images
            images = CatalogImage(catalog_id_s, band_type="MS", bbox=map(float, bbox.split(",")),proj=UTM_EPSG_code
                         ,pansharpen= False)
            imagew = CatalogImage(catalog_id_w, band_type="MS", bbox=map(float, bbox.split(",")),proj=UTM_EPSG_code
                         ,pansharpen= False)
        # calculate mean of coastal band to see if image is defective
            mean_coastals = images[1,:,:].read().mean()
            mean_coastalw = imagew[1,:,:].read().mean()

        # use image if not defective, otherwise
            if mean_coastals > 200:
                pass
            else:
                iS=iS+1

            if mean_coastalw > 200:
                pass
            else:
                iW=iW+1

            if mean_coastalw > 200 and mean_coastals > 200:
                imagestatus = 'fine'


        return selections,selectionw,images,imagew

# Calculate the difference in NDVI based on summer and winter image. Also returns a masked ndvi map for the summer image
def NDVIdiff(imagesummer,imagewinter,parkshape_utm):
    # get coordinates of the polygon so it can be plotted over the image.
    x,y = parkshape_utm.buffer(20).exterior.xy

    # resize polygon and plot polygon over image

    # subtract minimal values from utm polygon x and y to set 0,0 point as start
    x1 = np.subtract(x, min(x))
    y1 = np.subtract(y, min(y))

    # devide the x and y coordinate of the polygon by the size of the image to match both sizes
    x2 = np.divide(x1,max(x1)/imagesummer.shape[2])
    y2 = np.divide(y1,max(y1)/imagesummer.shape[1])

    poly = Polygon([[x2[i], y2[i]] for i in range(len(x2))])

    #summer first
    red = imagesummer[4,:,:].astype(np.float32)
    nir = imagesummer[6,:,:].astype(np.float32)
    ndvis = np.clip((nir - red)/(nir + red), -1, 1)

    #Winter next
    red = imagewinter[4,:,:].astype(np.float32)
    nir = imagewinter[6,:,:].astype(np.float32)
    ndviw = np.clip((nir - red)/(nir + red), -1, 1)

    #summer
    polygon = [(x2[i], y2[i]) for i in range(len(x2))]
    width = ndvis.shape[1]
    height = ndvis.shape[0]

    img = Image.new('L', (width, height), 0)
    ImageDraw.Draw(img).polygon(polygon, fill=1)
    mask = np.array(img)

    mask = mask == 1

    ndvis_test = copy.copy(ndvis)
    ndvis_test = np.flipud(ndvis_test)
    ndvis_test[mask == 0] = 0
    #filter out the water and roads
    ndvis_test[ndvis_test<0.1] = 0

    #winter
    ndviw_test = copy.copy(ndviw)
    ndviw_test = np.flipud(ndviw_test)
    ndviw_test[mask == 0] = 0
    #filter out the water and roads
    ndviw_test[ndviw_test<0.1] = 0

    #difference
    ndvi_diff = ndvis_test - ndviw_test

    # invert the mask so the outer boundaries are True
    mask = mask == 0
    # now calculate the mean NDVI values within the mask
    ndvims = np.ma.array(ndvis_test, mask=mask)
    ndvimw = np.ma.array(ndviw_test, mask=mask)
    ndvimd = np.ma.array(ndvi_diff, mask=mask)

    meansummer=ndvims.mean()
    meanwinter=ndvimw.mean()
    meandiff=ndvimd.mean()

    return meansummer,meanwinter, meandiff, ndvims
