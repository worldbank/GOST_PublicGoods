# Building Footprint Machine Learning

CREDITS: First scripts elaborated by Charles Fox and Sarah Antos for Bamako and then modified by Alex Chunet for easy and intuitive replication in Bujumbura and Harare. This analysis is currently being scaled up to other West African cities. An improved version of the methodology is currently being elaborated.

The objective of this project was to see whether we could train a machine learning algorithm to classify building footprints by income bracket using a small sample of training data.

The inputs were:
- A DigitalGlobe generated shapefile layer of all building footprints in Bamako
- A collection of simple shapefiles identifying low, medium, and high income residential regions

The analytical process can be summarized as follows:
1) Acquire DigitalGlobe/Maxar building footprint layer for AOI
2) Define training areas - informal areas / medium income / high income / commercial / industrial
3) Run Step 1 to build building footprint dataset with all necessary characteristics
4) Run Step 2 to train the model and make predictions

**Step 1 - Bamako Building Attribution** - This script attaches statistics about each footprint to
a GeoPandas GeoDataFrame describing each building footprint. It looks at various
characteristics, including:
- area of building footprint
- proximity and characteristics of 5 nearest neighbours
- proximity and characteristics of 25 nearest neighbours
- count of other buildings within 25 / 50 / 100 meters Depending on the size of the footprints file, it can take quite a while to run on a single processor. However, the script
has been designed to support multi-threading, so the more cores you have at your disposal, the faster it will run.

**Step 2 - Bamako Machine Learning** - Before running this script, we assume that you, the
user, have successfully run Step 1, and that you now have a building footprints layer
complete with statistics about each footprint. We also assume you generated some shapefiles
as training areas, which intersect some of the buildings and demarcate which income bucket
the buildings fall into.

This script takes the attributed building footprints layer, and the training area shapefiles, and trains 20 machine learning models concurrently using H2O's AutoML (see:
http://docs.h2o.ai/h2o/latest-stable/h2o-docs/automl.html). This technology has been selected mainly for ease of use - GOST is not specialized in machine learning deployment. The script will generate 20 models, select the best (in terms of accuracy of prediction - see aml.leaderboard) and then use this best model to generate predictions for the whole layer.

Users are encouraged to generate confusion matrices as a means of testing the accuracy of the process thereafter. Visualization of the final predictions can be easily done in QGIS. Note that, due to the limitations of H2O, it is not recommended to try to generate more than 100,000 predictions at a time; however, it is easy enough to break up large building footprint shapefiles into chunks, generate predictions for these chunks, and merge them afterwards, as demonstrated in the Step 2.

**Sample Output**
Eventually, if all goes well, the user should end up with predictions at the individual building
level, per the image below:

![Sample output](https://user-images.githubusercontent.com/35847289/175666406-35e77816-3464-4a60-9b74-619c2a527af4.png)
