
### home made functions
import nice_functions as nf

reload(nf)

import get_OSM_polygons as getOSM

reload(getOSM)

### other libraries

from PIL import Image, ImageDraw
import copy

import pandas as pd
import numpy as np

from shapely.ops import transform
from shapely.geometry import mapping, Polygon, box, shape
import matplotlib.pyplot as plt

import random

import fiona

import pickle

from gbdxtools import Interface
from gbdxtools.task import env
from gbdxtools import CatalogImage

from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

from sklearn.model_selection import train_test_split

from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier

from sklearn.tree import DecisionTreeClassifier

gbdx = Interface()


def get_classifier(city,buffer_size, min_size, sample_size):

    ################ this takes some time

    label_all = np.array([])
    data_all = np.array([])
    
    # labels of types to loop over searching for polygons
    labels = ['forest', 'grass', 'water' , 'building']

    # set minimal sizes for polygons of different types
    dict_size = {'forest': 1 , 'grass': .6, 'water': 1, 'building': .8} 

    # create empty dataframes that will be filled 
    geom_list_selection_all = []
    selection_all = pd.DataFrame()

    # loop over OSM types and get polygons 
    for label in labels: 

        selection, geojson_select, geom_list_selection,UTM_EPSG_code,project_utm,project_wgs = getOSM.get_OSM_polygons(city = city, type_query = label,min_size = dict_size[label])


        geom_list_selection_all.extend(geom_list_selection[1:sample_size])

        selection_all = selection_all.append(selection[1:sample_size])

    ################ this takes some time    

    selection_all = selection_all.reset_index().drop('index',axis = 1)
    
    # set classes for each polygon type
    dict_type = {'Forest': 1,'Wood': 1,'Nature Reserve': 3,'Wetland': 1, 'Grass': 2, 'Farmland': 2, 'Water': 3, 'Building': 4, 'Theatre': 0}

    # A list of "random" colors (for a nicer output)
    COLORS = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941"]
    


    # load multipolygon type from pickle for check
    multipolygon_type = pickle.load( open( "/home/gremlin/GGCW_tools_git/Pickle/multipolygon_type.p", "rb" ) )

    ################ this takes some time
    
    for objects in selection_all.index:

        print '-----------------------------------------\n'
        ### setting a buffer can remove the polygon or make it into a multipolygon, both are unusable so check if this is the case
        park_utm = transform(project_utm, geom_list_selection_all[objects])  # apply projection


        # perform check # get x y coordinates of polygon and set a buffer if polygon is large enough
        if (type(park_utm.buffer(buffer_size)) == multipolygon_type) | (park_utm.buffer(buffer_size).area == 0):

            message = "is multipolygon"


            print message + ' object: ' + str(objects) +'\n'

        else: 

            x,y = park_utm.buffer(buffer_size).exterior.xy

            park_buffer_wgs = transform(project_wgs,park_utm.buffer(buffer_size))  # apply projection  


            # get wgs projected x,y coordinates and create bounding box for image aquisition
            x_wgs,y_wgs = park_buffer_wgs.exterior.xy

            bbox_park_area_float = min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)

            bbox_park_area = str([min(x_wgs), min(y_wgs), max(x_wgs), max(y_wgs)])

            bbox_park_area_str = nf.listToStringWithoutBrackets(bbox_park_area)


            # convert bounding box to well known format usable by GBDX tools
            bbox_wkt = box(*bbox_park_area_float).wkt


            selection_images = nf.image_query_check(bbox_wkt,park_utm,buffer_size,multipolygon_type,project_wgs,x_wgs,y_wgs)

            if not selection_images.empty:
                # set park bounding box 
                bbox =  bbox_park_area_str

                # set catalog id from selection
                catalog_id =  selection_images.id[0]

                # collect image
                image = CatalogImage(catalog_id, band_type="MS", bbox=map(float, bbox.split(",")),
                                     proj=UTM_EPSG_code,pansharpen=False)


                #create array from GBDX image
                image_array = image[:,:,:].read()

                # get second band to see if image is defective (some images show only black)
                # use image if not defective, otherwise 
                if (image_array.size == 0):

                    print('no image')

                elif (image_array[1,:,:].min() != 0) :


                    print 'image for ' + selection_all.OSM_id[objects] + ' is good object: ' + str(objects )
                    print image.ipe_metadata["image"]["acquisitionDate"]
                    print image.cat_id
                    print selection_all.item_type[objects]


                    # resize polygon and plot polygon over image
                    # subtract minimal values from utm polygon x and y to set 0,0 point as start 
                    x1 = np.subtract(x, min(x))
                    y1 = np.subtract(y, min(y))

                    # devide the x and y coordinate of the polygon by the size of the image to match both sizes 
                    x2 = np.divide(x1,max(x1)/image.shape[2])
                    y2 = np.divide(y1,max(y1)/image.shape[1])


                    n_bands, rows, cols  = image.shape

                    # calculate total cells for each class by masking and setting pixel values to 1

                    # create sequence of edited x and y coordinates, widht and heigth  for use in ImageDraw function
                    polygon = [(x2[i], y2[i]) for i in range(len(x2))]
                    width = image.shape[2]
                    height = image.shape[1]

                    # convert polygon coordinates to raster/array values using ImageDraw
                    img = Image.new('L', (width, height), 0)
                    ImageDraw.Draw(img).polygon(polygon, fill=dict_type[selection_all.item_type[objects]])
                    # convert image to array and set as mask
                    mask = np.array(img)


                    # flip the array for matching with the mask array
                    image_array_flipped = np.fliplr(image_array[:,:,:])
                    reshaped_data = image_array_flipped.reshape(8,(rows*cols))
                    reshaped_label = mask.reshape(1,(rows*cols))

                    # check if this is the first iteration, if so add the first data set otherwise:
                    # append the new image data to the other data
                    if data_all.size == 0:

                        data_all = reshaped_data 

                    else: 

                        data_all = np.concatenate((data_all,reshaped_data), axis = 1)

                    label_all = np.append(label_all,reshaped_label)

                    # Two subplots, the axes array is 1-d
                    f, axarr = plt.subplots(1,2)
                    axarr[0].imshow(mask)
                    axarr[1].imshow(image_array_flipped[1])


                    plt.show()

    #                 plt.imshow(mask)
    #                 plt.show()

                    print label_all.shape
                    print data_all.shape


                else:

                    print 'image defective' 
                    # move to next without doing analysis

            else:

                print 'no image' 
                # move to next without doing analysis



    ################ this takes some time


    ### Remove pixels without class

    label_all_no0 = label_all[label_all != 0]
    data_all_no0 = data_all[:,label_all != 0]
    
    

    data = pd.DataFrame(data_all.T)


    data['class'] = label_all
    
    # Calculate NDVI

    # 8-band (0:Coastal, 1:Blue, 2:Green, 3:Yellow, 4:Red, 5:Red-edge, 6:NIR1, 7:NIR2) Multispectral

    # ndvi = (nir - red)/(nir + red)

    # EVI = 2.5 * ( nir - red ) / ( nir + 6.5 * red - 7.5 * blue+ 1.0 )

    data['ndvi'] = (data[6] - data[4])/(data[6] + data[4])

    data['EVI'] = 2.5 * (data[6] - data[4]) / (data[6] + 6.5 * data[4] - 7.5 * data[1] + 1 )


    data['water_index'] = (data[7] - data[0]) / (data[7] + data[0])


    X = data.iloc[:,0:8]

    y = data['class']


    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    # initialize search space (as a library!)
    param_grid = {
         'max_depth':[20,10],'n_estimators':[40]
    }

    gs = make_pipeline(     StandardScaler(), 
                            GridSearchCV(RandomForestClassifier(min_samples_leaf=2),
                            param_grid = param_grid,
                            cv = 2,
                            refit = True,
                            n_jobs = 1,
                            verbose = 0))

    # Instantiate random forest. You can specify default parameters here.
    # These parameters are not being optimized.


    # initialize grid search
    #gs = GridSearchCV(rf, param_grid, verbose=2)#,scoring='roc_auc')

    gs.fit(X_train,y_train)
    
    return gs