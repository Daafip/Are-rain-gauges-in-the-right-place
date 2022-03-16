import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from geopy.distance import geodesic

buildings = pd.DataFrame()

def distance_closest_house(station, buildings=buildings):
    '''
    Supply a data frame (df) in the format including atleast "station_latitude" and "station_longitude"
    Supply another including buildings with geometry.
    Returns closest station in area.
    
    '''
    # for every station create a list and dict
    lst = []
    # fany print update
    # loop over every other station
    for building in buildings['geometry']:
        # compute the distance
        distance = station.distance(building)*10**5
        # store that distance
        lst.append(distance)
    # now we've computed all the distances, take the minimum
#     d_min = min(lst)
    # if the closest station if further than the threshold
    return sorted(lst)[0]