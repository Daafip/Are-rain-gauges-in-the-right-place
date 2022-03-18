import numpy as np
import rasterio
import matplotlib.pyplot as plt

"""
Bulk of code to run 
"""

class Raster_Location_OnSlope:

    """\
        Main code to run all location on slope functions\
        Takes df: most importantly must have a raster ID\
        \
        """

    def __init__(self, df, lst_UK, n_px):
        self.df = df  # main data frame holding information
        self.lst_UK = lst_UK  # list containing the names of the rasters in the UK
        self.n_px = n_px
        self.label = {"ew":"East to west",
                    "ns": "North to south",
                    "nw-se": "North-west to south-east",
                    "ne-sw": "North-east to south-west"}


    def get_array_given_orientation(self, n, kind, clip=False): 
        """Takes index of gauges df, retrieves the file, return array of height data in the specified orientation"""
        raster_id_n = self.df.raster_id.iloc[n] 
        filepath = self.lst_UK[raster_id_n]

        with rasterio.open(f"Cropped Data/{self.lst_UK[raster_id_n]} - {n}.adf") as cropped:
            array = cropped.read(1)
            coords = rasterio.transform.rowcol(cropped.transform, self.df.station_lo.iloc[n], self.df.station_la.iloc[n])
        if kind == "ew":
            if clip:
                new_lst = []
                for i in array[coords[0],:]:
                    if i >= 0:
                        new_lst.append(i)
                    else:
                        new_lst.append(-1)
                return new_lst, coords[0]
            else:
                return array[coords[0],:], coords[0]
            
        elif kind == "ns":
            if clip:
                new_lst = []
                for i in array[:,coords[1]]:
                    if i >= 0:
                        new_lst.append(i)
                    else:
                        new_lst.append(-1)
                return new_lst, coords[1]

            else:
                return array[:,coords[1]], coords[1]
        elif kind == "nw-se" or kind == "ne-sw":
            if kind == "nw-se": 
                diagonal_array = np.diagonal(array)
            else:
                diagonal_array = np.fliplr(array).diagonal()
            if clip:
                new_lst = []
                for i in diagonal_array:
                    if i >= 0:
                        new_lst.append(i)
                    else:
                        new_lst.append(-1)
                return new_lst, coords[1]

            else:
                return diagonal_array, coords[1]

            
        else:
            ### no kind given
            print("incorrect kind")

    
    def distance_between_px(self, n):
        """need resolution between each pixe,  we know the resolution is 3 arc seconds, verified here https://www.opendem.info/arc2meters.html"""
        one_arc_second = np.cos(self.df.station_la.iloc[n] * np.pi / 180) * (1852/60)
        return one_arc_second * 3

    
    def slope_steepness(self, n, kind):
        """\
        Given an index, direction and length to consider, calculate: 
        rise / run * 100%
        rise = 50m per contour
        run  = computed distance
        Note in order of direction, meaning a slope ns with positive implies north is higher than south, 
        this also allows us to get the direction of slope. \
        """
        # get the array in a given orienteation
        array, location = self.get_array_given_orientation(n, kind)
        
        # crop down to find threshold:
        array_cropped = array[location-self.n_px:location+self.n_px+1]
        
        # remove any cut out pixels < -3000
        array_cropped = array_cropped[array_cropped >= -3000]

        return self.slope_steepness_from_array(n, array_cropped)

        
    def slope_steepness_from_array(self, n, array):
        """uses average method"""
        diff = []
        for index, item in enumerate(array):
            if index < len(array) - 1:
                rise = item - array[index+1]
                diff.append(rise)
        ### avoid divide by 0 error posed when doing steps of 1
        if (len(diff)-1) > 0:
            return ((sum(diff)/(len(diff)-1)) / self.distance_between_px(n)) * 100
        else:
            return (sum(diff) / self.distance_between_px(n)) * 100

    
    def get_array_step(self, n, kind, step_size, direction, rolling):
        """Given an index, direction and stepsize, returns the array"""
        # get the array in a given orienteation
        array, location = self.get_array_given_orientation(n, kind)

        self.array_length = len(array)//2
        # rolling means we dont average, make sure this is set correctly repending how you use the function
        if direction == "left":
            # go 'left' in the array
            if rolling:
                # rolling takes steps in the width of the stepsize, used for changes in slope steepness
                return array[location-self.n_px-step_size : location-self.n_px-step_size+2]
            else:
                # not rolling: used for averages slope location
                return array[location-self.n_px-step_size : location-self.n_px+1]

        elif direction == "right":
            # go 'right' in the array
            if rolling:
                return array[location+self.n_px+step_size-1 : location+self.n_px+step_size+1]
            else: 
                # not rolling: used for averages slope location
                return array[location+self.n_px : location+self.n_px+step_size+1]

        
        else:
            raise ValueError("Direction undefined")


    def get_location_on_slope(self, n):
        """For a given location determins where its situated.
        Looks at a spacing the same pixel width as initial steepness: main parameter therfore is self.n_px
        It looks left and right the same amount. This is an array width self.n_px shifted therefore use 'false'.
        Then runs it through if statements to look at the type."""

        lst = []
        for kind in self.label:
            initial_steepness = self.slope_steepness(n, kind)
            left = self.slope_steepness_from_array(n,self.get_array_step(n, kind, self.n_px, "left", False))
            right = self.slope_steepness_from_array(n,self.get_array_step(n, kind, self.n_px, "right", False))

            def determin_label(left, right):
                if left < 0 and right < 0: 
                    return "Slope"
                elif left > 0 and right > 0:
                    return "Slope"
                elif left < 0 and right > 0:
                    return "Ridge"
                elif left > 0 and right < 0: 
                    return "Valley"
                elif left == 0 == right:
                    return "Flat" 
                else:
                    return None
            
            location = determin_label(left, right)

            # print(location)

            # resurive to compare it to the steepness of the slope when constant without rewriting if statement
            if location is None:
                if left == 0:
                    location = determin_label(initial_steepness, right)
                elif right == 0:
                    location = determin_label(left, initial_steepness)
            if location is None:
                print(left,right)
            lst.append(location)
        return lst   

    def run_location_on_slope(self):
        """Takes list from get_location_on_slope and applies it to a df"""
        return self.df.level_0.apply(self.get_location_on_slope)
        


    """
    Look at this later, not in get_array_step() we went rolling for this part but not that
    """


    def change_slope_steepness(self, n, kind, step_size, threshold):
        """\
        Given index, direction step_size and threshold(in %): determins how when the slope steepness changes at intervals of the stepsize,\
            unitl the absolute change is more than the theshold. 
        \
        """
        # We want to make sure we get segements the sizes of step_size, not larger
        self.rolling = True
        # get the initial steepness
        initial_steepness = self.slope_steepness(n, kind)
        print(f'initial.....{initial_steepness}')

        # get first step: (we want rolling steps so therfore true)
        first_left, first_right = self.get_array_step(n, kind, step_size, "left", True), self.get_array_step(n, kind, step_size, "right", True)

        # remove any cut out pixels < -3000 & calculate steepness using function above
        left = self.slope_steepness_from_array(n, first_left[first_left >= -3000])
        right = self.slope_steepness_from_array(n, first_right[first_right >= -3000])
        
        # set first run
        count_left, count_right = 0, 0

        # first look left
        # (We check that the slope hasn't changed) and (we aren't outside the array)
        while abs(left - initial_steepness) < threshold and (count_left + self.n_px + step_size) < self.array_length:
            # update count to new position
            count_left += step_size
            print(f'left={left}, #{count_left}')
            # given the new position, update the 
            array =  self.get_array_step(n, kind, count_left, "left", True)
            # check list is not empty
            if len(array[array >= -3000]) == 0:
                # otherwise break
                break
            left = self.slope_steepness_from_array(n, array[array >= -3000])

        # then look right; same principle
        while abs(right - initial_steepness) < threshold and (count_right + self.n_px + step_size) < self.array_length:
            print(self.array_length)
            count_right += step_size
            print(f'right={right}, #{count_right}')
            array =  self.get_array_step(n, kind, count_right, "right", True)

            right = self.slope_steepness_from_array(n, array[array >= -3000])

        return count_left, count_right

        # # remove any cut out pixels < -3000
        # array_cropped = array_cropped[array_cropped >= -3000]
        
        # #average_method 
        # diff = []
        # for index, item in enumerate(array_cropped):
        #     if index < len(array_cropped) - 1:
        #         rise = item - array_cropped[index+1]
        #         diff.append(rise)





    
class Slope_Steepness:

    def __init__(self, df, lst_UK, n_px):
        self.df = df  # main data frame holding information
        self.lst_UK = lst_UK  # list containing the names of the rasters in the UK
        self.n_px = n_px
        self.label = {"ew":"East to west",
                    "ns": "North to south",
                    "nw-se": "North-west to south-east",
                    "ne-sw": "North-east to south-west"}

    def run_get_raster_pixel(self, n):
        """Takes index of gauges retrieves the file, finds location and returns height"""
        raster_id_n = self.df.raster_id.iloc[n] 
        filepath = self.lst_UK[raster_id_n]  
        
        with rasterio.open(f"Cropped Data/{self.lst_UK[raster_id_n]} - {n}.adf") as cropped:
            array = cropped.read(1)
            
            ### own boilerplate code
            # bounds = cropped.bounds
            # loc_lo = int(((-bounds[0] + self.df.station_lo.iloc[n]) * array.shape[0]) / (-bounds[0] + bounds[2]))
            # loc_la = int(((-bounds[1] + self.df.station_la.iloc[n]) * array.shape[1]) / (-bounds[1] + bounds[3]))
            
            # there is also a built in method... lets use that
            coords = rasterio.transform.rowcol(cropped.transform, self.df.station_lo.iloc[n], self.df.station_la.iloc[n])
            
            ### should be doable with the build in sample but i couldn't get it to work
    #         out = rasterio.sample.sample_gen(cropped,(self.df.station_lo.iloc[n],self.df.station_la.iloc[n]))
    #         out = cropped.sample((self.df.station_lo.iloc[n],self.df.station_la.iloc[n])
        return array[coords[0],coords[1]]

    def get_raster_pixel(self): 
        """applies the function run_get_raster_pixel to the df"""
        return self.df.level_0.apply(self.run_get_raster_pixel)

    def crop_raster(self, n):
        """Takes index of gauges, computes mask, cuts it out by copying and adjusting meta data and saving it"""
        raster_id_n = self.df.raster_id.iloc[n] 
        filepath = self.lst_UK[raster_id_n]
        ## open the wanted cell and load the data, creating a mask
        with rasterio.open(f"DEM data/{filepath}/{filepath}/w001001.adf") as src:
            crop_image, crop_transform = rasterio.mask.mask(dataset=src,
                                                            shapes=self.df.query(f"index == {n}").geometry,
                                                            crop=True)
            crop_meta = src.meta.copy()

        ## reassign new meta data 
        crop_meta.update({"driver": "GTiff",
                        "height": crop_image.shape[1],
                        "width": crop_image.shape[2],
                        "transform": crop_transform})

        ## save with new meta data
        with rasterio.open(f"Cropped Data/{self.lst_UK[raster_id_n]} - {n}.adf","w",**crop_meta) as dest:
            dest.write(crop_image)

    def slope_steepness_all_orientations(self, n): 
        """Function to check slope steepness in all directions"""
        lst = []
        for orientation in self.label:
            lst.append(self.slope_steepness(n,orientation))
        return lst

    def slope_steepness_max(self): 
        """Function to return the steepest slope of the given directions"""
        def run_slope_steepness_max(n):
            lst_abs = []
            lst = []
            for orientation in self.label:
                percentage = self.slope_steepness(n,orientation)
                lst_abs.append(abs(percentage))
                lst.append(percentage)
            return round(lst[lst_abs.index(max(lst_abs))],2)
        return self.df.level_0.apply(run_slope_steepness_max)


    def slope_steepness_max_abs(self): 
        """Function to return the absolute value of steepest slope of the given directions"""
        def run_slope_steepness_max_abs(n):
            lst_abs = []
            lst = []
            for orientation in self.label:
                percentage = self.slope_steepness(n, orientation)
                lst_abs.append(abs(percentage))
                lst.append(percentage)
            return round(max(lst_abs),2)
        return self.df.level_0.apply(run_slope_steepness_max_abs)


    def slope_steepness_max_direction(self): 
        """Function the direction in which the slope is steepest"""
        def run_slope_steepness_max_direction(n):
            label = ["ew","ns","nw-se","ne-sw"]
            lst_abs = []
            for orientation in label:
                percentage = self.slope_steepness(n,orientation)
                lst_abs.append(abs(percentage))
            return label[lst_abs.index(max(lst_abs))]    
        return self.df.level_0.apply(run_slope_steepness_max_direction)


    
class plotting:

    """\
    Main code to plot raster data\
    Inputs:
    - df: most importantly must have a raster ID\
    - lst_UK: list containing the names of rasterfiles e.g. "n50e000_dem"\
    - n_px: Number of pixels to be accounted for when plotting.
    \
    """

    def __init__(self, df, lst_UK, n_px):
        self.df = df  # main data frame holding information
        self.lst_UK = lst_UK  # list containing the names of the rasters in the UK
        self.n_px = n_px # amount of pixels to use
        self.label = {"ew":"East to west",
                    "ns": "North to south",
                    "nw-se": "North-west to south-east",
                    "ne-sw": "North-east to south-west"}


    def plot_cropped_raster(self, n, plot_hist=False, zoom_in=False, save=False):
        """Takes index of gauges DF retrieves the file and plots it."""
        raster_id_n = self.df.raster_id.iloc[n] 
        filepath = self.lst_UK[raster_id_n]
        
        # creating figures:    
        fig1, ax1 = plt.subplots(1, figsize=(10, 10))

        if plot_hist:
            fig2, ax2 = plt.subplots(1, figsize=(7, 5))  
        
        with rasterio.open(f"Cropped Data/{self.lst_UK[raster_id_n]} - {n}.adf") as cropped:
            plot = rasterio.plot.show(cropped, ax=ax1,cmap="terrain")
            if plot_hist:
                rasterio.plot.show_hist(cropped, ax=ax2, bins=10, 
                                        lw=0.0, stacked=False, alpha=0.9,
                                        histtype='stepfilled', title="Histogram")
        # post plotting
        fig1.colorbar(plot.get_images()[0], ax=ax1, shrink=0.75, label="Height(m)")
        ax1.plot(self.df.station_lo.iloc[n], self.df.station_la.iloc[n], "ro", markersize=5, label="Rain gauge")
        ax1.set_title(f'Surroundings of {self.df.station_fi.iloc[n]} at ({self.df.station_lo.iloc[n]},{self.df.station_la.iloc[n]})',
                    pad=18);
        
        if zoom_in:
            # probably a better way to dit it but this works, as only 3 decimal palces this gives an esimate of accuracy
            ax1.vlines([self.df.station_lo.iloc[n]+ 0.0005, self.df.station_lo.iloc[n]-0.0005],
                    self.df.station_la.iloc[n] - 0.05,self.df.station_la.iloc[n] + 0.05, 
                    color="red",linestyles="dashed")
            ax1.hlines([self.df.station_la.iloc[n]+ 0.0005, self.df.station_la.iloc[n]-0.0005],
                    self.df.station_lo.iloc[n] - 0.05,self.df.station_lo.iloc[n] + 0.05, 
                    color="red",label="accuracy limit",linestyles="dashed")
            ax1.set_xlim([self.df.station_lo.iloc[n]+ 0.005, self.df.station_lo.iloc[n]-0.005])
            ax1.set_ylim([self.df.station_la.iloc[n]+ 0.005, self.df.station_la.iloc[n]-0.005])
            
        if plot_hist:
            ax2.legend("")
            ax2.set_xlabel("Height above MSL")
        
        ax1.legend(bbox_to_anchor=(1.25, +0.04))
        if save:
            ax1.figure.savefig(f"Plots/{n} - Surroundings of {self.df.station_fi.iloc[n]}.jpg")



    def plot_gauge(self, n, kind="ew", clip=True, save=False):
        """Function which takes index of a gauge and what kind of plot wanted"""    
        instance = Raster_Location_OnSlope(self.df,self.lst_UK,self.n_px)
        if kind in self.label:
            array, gauge_loc = instance.get_array_given_orientation(n, kind, clip)
            # plotting
            fig, ax = plt.subplots(1, figsize=(6, 6))
            ax.plot(array, label=self.label[kind])
            ax.axhline(-1,color='g',label="End map", alpha=0.5)
            ax.plot(gauge_loc,Slope_Steepness(self.df,self.lst_UK,self.n_px).run_get_raster_pixel(n),"ro",label="Gauge")
            ax.legend(bbox_to_anchor=(1.00, -0.14))
            title = f"Profile of slope around gauge #{n}"
            ax.set_title(title)
            ax.set_xlabel("Pixels")
            ax.set_ylabel("Height[m]")
            if save:
                fig.set_figheight(7)
                fig.set_figwidth(7)
                ax.legend()
                ax.figure.savefig(f'Figure/{title}{kind}.jpg')
                
        # recursively calls all the options & plots
        elif kind == "all":
            # plotting fancy
            kind_array = [["ew","ns"],["nw-se","ne-sw"]]
            fig, ax = plt.subplots(2,2,constrained_layout=True, figsize=(5, 5))
            title = f"Profile of slope around gauge #{n}"
            fig.suptitle(title,y=1.04,fontsize=16)
            # fig.tight_layout()
            for i in range(2):
                for j in range(2):
                    kind = kind_array[i][j]
                    array, gauge_loc = instance.get_array_given_orientation(n, kind, clip)
                    ax[i,j].plot(array)
                    ax[i,j].axhline(-1,color='g',label="End map", alpha=0.5)
                    ax[i,j].plot(gauge_loc,Slope_Steepness(self.df,self.lst_UK,self.n_px).run_get_raster_pixel(n),"ro",label="Gauge")
                    ax[i,j].set_title(f"Slope {self.label[kind]}")
                    ax[1,j].set_xlabel("Pixels")
                    ax[i,0].set_ylabel("Height[m]")
                    ax[0,j].get_xaxis().set_visible(False)
                    ax[i,1].get_yaxis().set_visible(False)

            ax[1,1].legend(bbox_to_anchor=(1.00, -0.05))
            
            if save:
                plt.savefig(f'Figure/{title} - all.jpg');   
   
        else:
            print("incorrect kind")   


    def slope_steepness_all_orientations_varypxs(self, n, n_px): 
        """Function to check slope steepness in all directions"""
        lst = []
        for orientation in self.label:
            lst.append(self.slope_steepness_varypxs(n, n_px, orientation))
        return lst

    def slope_steepness_varypxs(self, n, n_px, kind):
        """ copy of prev but edited"""
        instance = Raster_Location_OnSlope(self.df,self.lst_UK,self.n_px)
        # get the araay
        array, location = instance.get_array_given_orientation(n, kind)
        # crop down to the length to consider
        array_cropped = array[location-n_px:location+n_px+1]
        
        # remove any cut out pixels < -3000
        array_cropped = array_cropped[array_cropped >= -3000]
        
        #average_method 
        diff = []
        for index, item in enumerate(array_cropped):
            if index < len(array_cropped) - 1:
                rise = item - array_cropped[index+1]
                diff.append(rise)
                
        return ( (sum(diff)/(len(diff)-1)) / instance.distance_between_px(n)) * 100


    def plot_varying_n_px(self, lst, subplots=True,save=False):
        """\
        Plots for a list of supplied stations, the resulst of slope steepness.
        Subplots: define seperate plots or not
        save: save or not 
        """
        for i in lst:
            lst = []
            for j in range(1,30):
                n_px = j

                ## again little jank but will do
                
                lst.append(self.slope_steepness_all_orientations_varypxs(i, n_px))
            if subplots:
                fig, ax = plt.subplots(1,figsize=(4,4))
                ax.plot(lst)
                legend = ["East to west","North to south","North-west to south-east","North-east to south-west"]
                ax.legend(legend,bbox_to_anchor=(1.8,0.35))
                ax.set_xlabel("pixels taken into account")
                ax.set_ylabel("Corresponding slope percentage")
                title = f"Varying pixels for gauge #{i}"
                ax.set_title(title)
                if save:
                    fig.set_figheight(7)
                    fig.set_figwidth(7)
                    ax.legend(legend)
                    ax.figure.savefig(f'Figure/{title}.jpg',bbox_extra_artists=["legend"])
                    
            else:
                plt.plot(lst)

        if not subplots:
        # set art for all 
            plt.legend(["East to west","North to south","North-west to south-east","North-east to south-west"],
                            bbox_to_anchor=(1,-0.15))
            plt.xlabel("pixels taken into account")
            plt.ylabel("Corresponding slope percentage")
            title = f"Overview of slope percentage varying pixel amounts"
            plt.title(title)      
            if save:
                plt.savefig(f"Figure/{title}.jpg")
    
    
    def get_plot(self, n):
        """Simple array plot"""
        filepath = self.lst_UK[n]
        with rasterio.open(f"DEM data/{filepath}/{filepath}/w001001.adf") as src:
            out = src.read(1)

        return out.clip(-200)
