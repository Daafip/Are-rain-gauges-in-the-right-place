import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
import geopandas
import time
from geopy.distance import geodesic

def distance_too_close(df, threshold):
    '''
    Supply a data frame (df) in the format including atleast "station_latitude" and "station_longitude"
    Returns a list of stations closer than the thereshold
    
    '''
    ## adds a location column to the dataframe
    df['Location'] = df[['station_longitude', 'station_latitude']].apply(tuple, axis=1)
    # empty list to store
    lst = []
    # loop over all the locations
    for loc1 in df['Location']:
        # fancy update print statement
        print(loc1,end="\r")
        # loop over all the locations again
        for loc2 in df['Location']:
            # no need to caculate the distance of 0
            if loc1 == loc2:
                pass 
            else:
                # calculate the distance beteween points
                distance = geodesic(loc1,loc2).km
                # if intrested we store it
                if distance < threshold:
                    # first check if the other point is in the list to avoid doubles
                    if (loc2,loc1,distance) in lst:
                        pass
                    else:
                        lst.append((loc1,loc2,distance))
                        # if not append
    # finally return a list with tuples containing coordinates and distances
    return lst

def tuple_to_csv(list_tuples, filename="tuple_to_csv-file"):
    """
    Takes a list of tuples in the form (station1, station2, distance apart) and writes the locations to a csv file
    Only input name, .csv is assumed
    This gives a format which is nice for GIS importing
    """
    lst_out = "station_latitude,station_longitude,distance\n"
    for i in list_tuples:
        lst_out += str(i[0][1]) + "," + str(i[0][0]) + "," + str(i[2]) + "\n" \
        + str(i[1][1]) + ","+ str(i[1][0]) + "," + str(i[2]) +"\n"
    with open('Output/' + filename + '.csv', 'w') as f:
        f.write(lst_out)
        
def distance_closest(df, threshold):
    '''
    Supply a data frame (df) in the format including atleast "station_latitude" and "station_longitude"
    Returns a list of stations where the closest if above the thereshold
    
    '''
    # creates column with tupel (lon,lat)
    df['Location'] = df[['station_longitude', 'station_latitude']].apply(tuple, axis=1)
    # list to store results
    lst_total = []
    # loop over all stations
    for loc1 in df['Location']:
        # for every station create a list and dict
        lst = []
        lookup_loc2 = {}
        # fany print update
        print(loc1,end="\r")
        # loop over every other station
        for loc2 in df['Location']:
            # skip itself
            if loc1 == loc2:
                pass
            else:
                # compute the distance
                distance = geodesic(loc1,loc2).km
                # store that distance
                lst.append(distance)
                # store the corresponding tuple of the station in de dict
                lookup_loc2[str(distance)] = loc2
        # now we've computed all the distances, take the minimum
        d_min = sorted(lst)[0]
        # if the closest station if further than the threshold
        if d_min > threshold:
            # and we haven't already stored the pair
            if (lookup_loc2[str(d_min)],loc1,d_min) in lst_total:
                pass
            else:
                # then append the tuple in the same form as before
                lst_total.append((loc1,lookup_loc2[str(d_min)],d_min))
    return lst_total
 

def run_buffer_spacing(df, name):
    """takes df, computes buffer around it, and computes area"""
    lst = []
    buffer = []
    for i in np.linspace(0.0005,0.5,num=50):
        instance = df.copy()
        warnings.filterwarnings("ignore")
        instance.geometry = instance.geometry.buffer(i)
        lst.append(instance)
        buffer.append(i)

    start = time.time()
    nmax = len(lst) # total amount of variations of buffering
    imax = len(lst[0]) # amount of stations
    out = np.zeros((nmax,imax))
    for n in range(nmax):
        print(f'n={n}',end="\n")
        for i in range(imax):
            print(i,end="\r")
            out[n,i] = sum(lst[n].query(f'index == {i}').overlay(lst[n].query(f'index != {i}'), how='intersection').geometry.area)
    end = time.time()
    print((end-start)/60)
    
    np.savetxt(f"Output/intersection-array-{name}.csv",out,delimiter=",")
    np.savetxt(f"Output/buffer-array-{name}.csv",buffer,delimiter=",")
    return out, lst, buffer
    
