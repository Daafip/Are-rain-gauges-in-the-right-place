# Are rain gauges in the right place?
Collection of code written as part of Bachelor Thesis Civil Engineering:  'Are raingauges in the right place?' - March 2022

Folder each consist of a different part of the research (in order of refactoring):


#### _Used in thesis:_
- Distance stations: loads in gauge locations, looks at the distances between rain gauge stations (basic intro)
- Difference Contour and DEM: Loads in results from contour and DEM analysis, merges the dataframes and creates plots (plots)
- Contours: Takes shapefile contour data in the format supplied by the  [Ordnance Survey](https://osdatahub.os.uk/) at 50m resolution and turns it into slope steepness data for specific points. In this case raingauge points but could be used at any given location in theory. 


#### _unused:_
- Distance building: Loads in shapedata for buildings and computes the distance to buidlings (exploring use of shapely)
