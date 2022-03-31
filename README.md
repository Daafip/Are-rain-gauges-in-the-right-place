# Are rain gauges in the right place?
Collection of code written as part of Bachelor Thesis Civil Engineering:  'Are raingauges in the right place?' - March 2022



#### _Folder each consist of a different part of the research (in alphabetical order)_
- Difference Contour and DEM: Loads in results from contour and DEM analysis, merges the dataframes and creates plots (plots).
- Distance building: Loads in shapedata for buildings and computes the distance to buidlings ( _Exploring use of shapely, unused in final report_).
- Distance stations: loads in gauge locations, looks at the distances between rain gauge stations (basic intro).
- Location on slope from DEM: takes the data computed in slopes on DEM and further processes this four categories of possible locations on slopes.
- Slope steepness from DEM: Takes void filled elevation (raster) data in format supplied by [hydroSHEDS](https://hydrosheds.org) at 3 arc-second resolution and turns it into slope steepness data.
- lope steepness from Contours: Takes shapefile contour data in the format supplied by the  [Ordnance Survey](https://osdatahub.os.uk/) at 50m resolution and turns it into slope steepness data for specific points. In this case raingauge points but could be used at any given location in theory. 
- 

#### _Data usage:_

All data used is open source and referenced where used. Hosted here to allow for easy use and compatibility when using the code. 

#### Packages:
environment.yaml contains all necessary packages to run the data, can be loaded into [Anaconda](https://www.anaconda.com/) to set up for correct use.
