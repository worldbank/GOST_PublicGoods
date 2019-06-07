### GOSTNets
This package contains a library and example notebooks for running network analysis leveraging a number of existing network analysis libraries: peartree, osmnx, and networkx

The package contains four folders:
1. GOSTNets - the library containng the functions for running the network analysis
2. SampleData - data for running examples found in Notebooks
3. Notebooks - sample notebooks for runnning analyses
4. Tutorial - introductory tutorial for getting to grips with the basics of the library. 

Every function contains a docstring which can be brought up in use to check the inputs for various functions. For example: 

```python
gn.edge_gdf_from_graph?
```

returns: 

```
Signature: gn.edge_gdf_from_graph(G, crs={'init': 'epsg:4326'}, attr_list=None, geometry_tag='geometry', xCol='x', yCol='y')
#### Function for generating a GeoDataFrame from a networkx Graph object ###
 REQUIRED: a graph object G
 OPTIONAL: crs - projection of format {'init' :'epsg:4326'}. Defaults to
           WGS84. Note: here we are defining the crs of the input geometry -
           we do NOT reproject to this crs. To reproject, consider using
           geopandas' to_crs method on the returned gdf.
           attr_list: list of the keys which you want to be moved over to
           the GeoDataFrame.
           geometry_tag - the key in the data dictionary for each edge which
           contains the geometry info.
           xCol - if no geometry is present in the edge data dictionary, the
           function will try to construct a straight line between the start
           and end nodes, if geometry information is present in their data
           dictionaries.  Pass the Longitude info as 'xCol'.
           yCol - likewise, determining the Latitude tag for the node's data
           dictionary allows us to make a straight line geometry where an
           actual geometry is missing.
 RETURNS: a GeoDataFrame object of the edges in the graph
#-------------------------------------------------------------------------#
```

These docstrings have been written for every function, and should help new and old users alike with the options and syntax.
