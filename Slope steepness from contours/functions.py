import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import warnings

##### Some sections can still be refactored to classes


"""
Loads in all the files:
"""


def bulk_import_shapefiles(file_prefix="ny",n=2,file_type="point"):
    """
    Allows  easy import of shapefiles, use the format as shown in the function.
    Returns two lists: one lst of the geopandas dataframe, other the centres of each grid.
    """
    lst_files = []
    lst_grid = []
    for i in range(n+1):
        # add prefix 0 to single digits
#         print(i,end="\r")
        if i < 10:
            index = f"0{i}"
        else:
            index = i
        try:
            # define path name
            filename = f"Contour Data/Contours/{file_prefix}/{file_prefix.upper()}{index}_{file_type}.shp"
            with open(filename):
                pass
            # read file
            shapefile = gpd.read_file(filename)
            shapefile.to_crs(epsg=4326,inplace=True)
            if file_type == "point":
                # indexing code to later find which shapefile we want per gauge: this is with centroids
                warnings.simplefilter("ignore")
                lst_grid.append(shapefile.dissolve().centroid[0])
            ## Inefficient Code but interesting none the less
#             if file_type is None:
#                 # indexing code to later find which shapefile we want per gauge: this is with the bounds of the contours
#                 boundaries = shapefile.dissolve().bounds.to_numpy()[0]
#                 lst_grid.append(Point((boundaries[0] + boundaries[2])/2,(boundaries[1] + boundaries[3])/2))
                
            
            # append the df to a list
            lst_files.append(shapefile)
            
        # NY31 doesn't exsist so we do it this way to deal with it
        except FileNotFoundError:
            lst_files.append(None)
            lst_grid.append(None)
#             print(f"Missing:{filename}")
            pass
    # finally return it 
    return lst_files,lst_grid


"""
This section finds and matches the closest grid. 
The following code is not the most ellegant way to do it, but it works. 
"""


def select_closest_grid(df, lst_grid_point):
    """
    Takes as df and a list of points to compare it to, finds grid refference closest to the gauge and returns that index.
    """
    def closest_grid_id(gauge):
        lst_distances = []
        # walk over all points
        for i in lst_grid_point:
            # check for nones
            if i is None:
                lst_distances.append(np.inf)
            else:
                # caculate distance between point and the gauge, appen to the list
                lst_distances.append(i.distance(gauge))
        # once all is done take the lowest
        closest = min(lst_distances)
        # return the index of the closest grid
        return lst_distances.index(closest)

    return df.geometry.apply(closest_grid_id)


def select_closest_dist(df, lst_grid_point):
    """
    Similar to select_closest_grid() but this time returns the distance, could be more efficient way but this works. 
    Takes as df and a list of points to compare it to, finds grid refference closest to the gauge and returns that distance.
    """
    def closest_gauge_dist(gauge):
        lst_distances = []
        for i in lst_grid_point:
            if i is None:
                lst_distances.append(np.inf)
            else:
                lst_distances.append(i.distance(gauge))
        closest = min(lst_distances)*10**5
        # same as before but returns the distance, not the index
        return closest

    return df.geometry.apply(closest_gauge_dist)

"""
"""
"""
This section computes the slope based on closest points on the contour
"""
"""
"""


class slope_steepness:
    
    """"All calculations to do with slope steepness""" 
    
    def __init__(self, df, lst_lines, lst_points, selected_n):
        self.df = df
        self.lst_lines = lst_lines
        self.lst_points = lst_points
        self.selected_n = selected_n
    
    
    """
    """
    """
    First done by distance to gauge
    """
    """
    """
    
    
    def by_distance_to_gauge(self):
        """
        Sorted by the distance from gauge to contours
        rise / run * 100%
        rise = 50m per contour
        run  = computed distance
        Takes averages and then turn into slope
        """
        def sorted_slope_steepness(index):
            d = self.sorted_distance_between_points_on_contours(index)
            return round((50 * 100) /(sum(d)/len(d)),2)

        return self.df.level_0.apply(sorted_slope_steepness)
    
    """Method 2 does it a little different, first computes the individual slopes then averages it, seems to work too but method 1' results seem better."""
    
    # def by_distance_to_gauge_method2(self):
    #     """
    #     Turns into slope and then averages
    #     """
    #     def sorted_slope_steepness(index):
    #         d = self.sorted_distance_between_points_on_contours(index)
    #         lst = []
    #         for i in d:
    #             if i > 0:
    #                 lst.append(50 / i * 100)
    #         return sum(lst)/len(lst)
    #     return self.df.level_0.apply(sorted_slope_steepness)
    
    
    """
    """
    """
    Then resorted by using the distance between contours
    """
    """
    """
    
    
    def by_closest_inter_distance(self):
        """
        Resorts all the points to 
        Takes averages and then turn into slope
        """
        def resorted_slope_steepness(index):
            d = self.resorted_points_on_contours(index)
            return round((50 * 100) /(sum(d)/len(d)),2)
        
        return self.df.level_0.apply(resorted_slope_steepness)

    """Method 2, see description above"""
    
    # def by_closest_inter_distance_method2(self):
    #     """
    #     Turns into slope and then averages
    #     """
    #     def resorted_slope_steepness(index):
    #         d = self.resorted_points_on_contours(index)
    #         lst = []
    #         for i in d:
    #             if i > 0:
    #                 lst.append(50 / i * 100)
    #         return sum(lst)/len(lst)
    #     return self.df.level_0.apply(resorted_slope_steepness)
    
    
    """
    """
    """
    All dependancy code:
    """
    """
    """


    def get_df_per_station(self, index):
        """
        Insert index of the gauges df and returns spotheights and contour dfs
        """
        station = self.df.iloc[index]
        spot_heights = self.lst_points[station.grid_id]
        contour_lines = self.lst_lines[station.grid_id]
        return spot_heights, contour_lines


    def select_closest_contour(self, index):
        """
        Takes gauge_id and returns df sorted by proximity of their centroids
        """
        ### station needs to be defined
        station = self.df.iloc[index]
        def distance_to_gauge(df):
            """Calculates distance to given station and adds to df column"""
            return df.distance(station.geometry)

        spot_heights, contourlines = self.get_df_per_station(index)
        contourlines['distance_to_gauge'] = contourlines.geometry.apply(distance_to_gauge)*10**5
        return contourlines.sort_values('distance_to_gauge')


    def distance_to_n_contours(self, index):
        """
        Repeats select_closest_contour() but for index of df, returns new df
        """
        contours = self.select_closest_contour(index).head(self.selected_n)
        point = self.df.iloc[index].geometry
        distance_contour = []
        for i in range(self.selected_n):
            line = contours.iloc[i].geometry
            distance_contour.append(point.distance(line.interpolate(line.project(point)))*10**5)
        return contours


    def closest_point_on_contour(self, index):
            """
            Takes the amount of contours specified by selected_n and finds the point closest to the gauge
            """
            contours = self.select_closest_contour(index).head(self.selected_n)
            rain_gauge = self.df.iloc[index].geometry

            def point_on_contour(line):
                """finds point on contour used to project and  adds to df column"""
                return line.interpolate(line.project(rain_gauge))
            
            contours["point_on_contour"] = contours.geometry.apply(point_on_contour)
            return contours


    def sorted_distance_between_points_on_contours(self, index):
        """
        Takes an index, returns the amount of points specified by self.selected_n
        """
        contour = self.closest_point_on_contour(index)
        points = contour.point_on_contour.to_numpy()
        distance = []
        
        for i in range(len(points)):
            if i == len(points) - 1:
                pass
            else:
                distance.append(round(points[i].distance(points[i+1])*10**5))
                
        return distance
    
    def resorted_points_on_contours(self, index):
        """
        Resorts all the points on the contours based on the distance between them
        """
        contour = self.closest_point_on_contour(index)
        points = list(contour.point_on_contour.to_numpy())

        # start with the furthest Node
        current_node = contour.point_on_contour.iloc[self.selected_n-1]

        # add this to the list of current nodes
        index_chosen_points = [self.selected_n-1]
        chosen_points = [current_node]
        chosen_distances = []

        # and remove the node from out list, replacing it with none to preserve indexing
        points[self.selected_n-1] = None

        #run untill 1 ellement left 
        while len(index_chosen_points) < len(points):

            # create a new list to store distances
            lst_distance = []

            # loop over all points
            for i in points:
                # check its not already chosen
                if i is not None:
                    # caculate distance from current to all others
                    lst_distance.append(round(current_node.distance(i)*10**5))
                else:
                    # again preserving indexes
                    lst_distance.append(np.inf)

            # take the minimum of the caculated nodes and find the index
            closest = min(lst_distance)
            # print(lst_distance)
            new_index = lst_distance.index(closest)

            # reassign new node and remove the current_node from the points list
            current_node = contour.point_on_contour.iloc[self.selected_n-1]
            points[new_index] = None

            # add this index to the selected list
            index_chosen_points.append(new_index)
            chosen_points.append(current_node)
            chosen_distances.append(round(closest))

        return chosen_distances
    
    
    """
    """
    """
    Plotting:
    """
    """
    """
    
    
    def plot_station(self, index, annotate=False):
        """
        Plots just a station with surroundings
        """
        station = self.df.iloc[index]
        base = self.lst_points[station.grid_id].plot(label= "Spot Heights",color="blue", alpha = 0.6)
        plt.plot(station.station_lo,station.station_la,"ro",label=f"{station.station_na} weather station")
        self.lst_lines[station.grid_id].plot(ax=base,alpha =0.1)
        plt.title(f"Distance to gridcell center:{station.grid_id_distance/1000:.3f}km ")
        plt.legend(bbox_to_anchor=(1, -0.05))
        
        if annotate:
            height = lst_points[station.grid_id].PROP_VALUE.to_numpy()
            x = lst_points[station.grid_id].geometry.x.to_numpy()
            y = lst_points[station.grid_id].geometry.y.to_numpy()
            for i, txt in enumerate(height):
                base.annotate(txt,(x[i],y[i]),label="Height above mean sea-level")

                
    def plot_closest_point_contours(self, index, annotate=False, spot_heights=False, distances=False):
        """
        Plots the given data till now
        """

        #calls main df and extracts needed data
        station = self.df.iloc[index]
        base = self.select_closest_contour(index).head(self.selected_n).plot(alpha =0.5)
        plt.plot(station.station_lo,station.station_la,"ro",label=f"{station.station_na} weather station")

        # adding the distance as title
        if distances:
            line = self.select_closest_contour(index).head(self.selected_n).iloc[0].geometry
            point = self.df.iloc[index].geometry
            distance_contour = point.distance(line.interpolate(line.project(point)))*10**5
            plt.title(f"Distance to closest contour: {distance_contour:.3f}m ")


        # plots spot heigts as wel as 
        if spot_heights:
            self.lst_points[station.grid_id].plot(ax=base, label= "contour_points",color="blue")

        # adds centroids and points on the contours
        if annotate:
            height = self.lst_points[station.grid_id].PROP_VALUE.to_numpy()
            x = self.lst_points[station.grid_id].geometry.x.to_numpy()
            y = self.lst_points[station.grid_id].geometry.y.to_numpy()
            for i, txt in enumerate(height):
                base.annotate(txt,(x[i],y[i]))
            height = self.lst_lines[station.grid_id].PROP_VALUE.to_numpy()

        # adds centroids and points on the contours    
        if self.selected_n < 31:
            self.select_closest_contour(index).head(self.selected_n).centroid.plot(ax=base,label="centroid",color="green")
            self.closest_point_on_contour(index).point_on_contour.plot(ax=base,label="Closest points on contours",color="purple")
        plt.legend(bbox_to_anchor=(1, -0.05))

        # returns for more plotting after
        return base
        
        
        """"""""" The end """""""""

    
        