import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
import math
import copy
from shapely.ops import transform
from shapely.geometry import mapping, Polygon, box, shape
from skimage.transform import resize
from skimage import io, filters, measure, color, morphology
from scipy import ndimage
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

from collections import namedtuple

from gbdxtools import Interface
from gbdxtools.task import env
from gbdxtools import CatalogImage

gbdx = Interface()

def Waterextract(image,parkshape_utm,cellsize):
    
#     Calculate the waterindex first
    coastal = image[0,:,:].astype(np.float32)
    nir2 = image[7,:,:].astype(np.float32)
    water_index = (nir2 - coastal)/(nir2 + coastal)
    # First, clean up any nan values
    water_index[np.isnan(water_index)] = 0

# Next, run a gaussian smoothing kernel over the image. This will smooth out localized noise in the water index
# by use a focal moving window.
    water_index_smoothed = filters.gaussian(water_index, preserve_range=True)

    # set the threshold
    binary_threshold = -0.65 

# Apply the threshold
    binary_water = water_index_smoothed >= binary_threshold
# Specify the minimum feature size in square meters, and then use info from the image metadata 
# to translate to grid cell count.

# Set the minimum feature size to 40 sq m
    min_feature_size_m2 = 40.

# From the image metadata, we can determine the area of a single grid cell
    cell_height_m = cellsize
    cell_area_m2 = cell_height_m**2

# Finally, use the cell size to convert the minimum feature size to grid cells
    min_feature_size_cells = np.round((min_feature_size_m2/cell_area_m2)).astype('int64')

# First, remove the small holes
    water_cleaned = morphology.remove_small_holes(binary_water, min_feature_size_cells)

# Then remove the small objects
    water_cleaned = morphology.remove_small_objects(water_cleaned, min_feature_size_cells, connectivity=2)
    
# create sequence of edited x and y coordinates, widht and heigth  for use in ImageDraw function

    #  first take the shapely file and write it to x an y.
    x,y = parkshape_utm.buffer(20).exterior.xy
    # subtract minimal values from utm polygon x and y to set 0,0 point as start 
    x1 = np.subtract(x, min(x))
    y1 = np.subtract(y, min(y))

# devide the x and y coordinate of the polygon by the size of the image to match both sizes 
    x2 = np.divide(x1,max(x1)/image.shape[2])
    y2 = np.divide(y1,max(y1)/image.shape[1])
    
    polygon = [(x2[i], y2[i]) for i in range(len(x2))]
    width = water_cleaned.shape[1]
    height = water_cleaned.shape[0]


# convert polygon coordinates to raster/array values using ImageDraw
    img = Image.new('L', (width, height), 0)
    ImageDraw.Draw(img).polygon(polygon, fill=1)
    # convert image to array and set as mask
    mask = np.array(img)

    # create boolean mask
    mask = mask == 1

    # copy layers for classification
    Water = water_cleaned

    # flip image and EVI to allign with polygon based mask
    Water = np.flipud(Water)


    # apply mask to test image
    Water = Water+1
    Water[Water==2] = 0
    Water[mask == 0] = 0
    
    return Water

# use Scipy to segment water poygons and measure some properties
def Watersegment(Water,cellsize):
    # Convert into Boolean
    Water_mask = Water >= 1
    # Use Scipy image to segment each polygon
    labeled_mask, num_labels = ndimage.label(Water_mask, structure=[[1,1,1],[1,1,1],[1,1,1]])
    clusters = measure.regionprops(labeled_mask, Water)
    
    # Create dataframe for water body properties
    Water_df = pd.DataFrame(columns=['id','Area','Eccentricity','Maj_Axis_Length','Min_Axis_Length','Perimeter'])
                                 
    for i in clusters:
        label=i.label
        ecc= i.eccentricity
        area= i.area*cellsize
        maja= i.major_axis_length*cellsize
        mina= i.minor_axis_length*cellsize
        per= i.perimeter*cellsize
   
    
   # load all metadata to dataframe
        Water_df = Water_df.append({'id':label ,'Area':area, 'Eccentricity':ecc,'Maj_Axis_Length':maja,'Min_Axis_Length':mina,'Perimeter':per},ignore_index=True) 
    return Water_df


def Wateredge(Water,cellsize, buffer):
    # Lets use a sobel filter to find the water edges
    sx = ndimage.sobel(Water, axis=0, mode='constant')
    sy = ndimage.sobel(Water, axis=1, mode='constant')
    sob = np.hypot(sx, sy)
    borders = sob >=1
    
# Use 2D convolution to create a buffer around the water edges
# first calculate number of cells needed for buffer based on cell size and desired buffer (20m)
    buf= int(buffer/cellsize)
    kernel = np.ones((buf, buf))
    
    result = np.int64(convolve2d(borders, kernel, mode='same') > 0)
    # Delete the water part so only the riperian zone is left
    riparian=result-Water
    
    return riparian

def Riparianveg(riparian,ndvi):
    # convert image to array and set as mask
    mask = np.array(riparian)

    
# create boolean mask
    mask = mask == 1

# copy layers for classification
    ndvi2 = copy.copy(ndvi)


# apply mask to test image
    ndvi2[ndvi2 >= 0.3] = 2
    ndvi2[ndvi2 < 0.3] = 1
    ndvi2[mask == 0] = 0
    maskcount = 1*mask
    
    vegetation = sum(sum(sum([ndvi2 == 2])))
    total = sum(sum(maskcount))

    percvegcover = np.divide(vegetation*100,total)
    
    return percvegcover