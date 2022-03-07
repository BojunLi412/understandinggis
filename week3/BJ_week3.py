#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 21:34:32 2021

@author: ulrica
"""



from math import sqrt
from geopandas import read_file
from rtree import index

def distance(x1, y1, x2, y2):   
    """
    * Use Pythagoras' theorem to measure the distance
    * This is acceptable in this application because the distances are so short
    *  and computational efficiency is a key consideration
    """
    return sqrt((x1-x2)**2+(y1-y2)**2)

pop_points = read_file("../data/gulu/pop_points.shp")
water_points = read_file("../data/gulu/water_points.shp").to_crs(pop_points.crs)
gulu_district = read_file("../data/gulu/district.shp").to_crs(pop_points.crs)

idx = index.Index()

# loop through each row and construct spatial index on the larger dataset (population)
for id, well in water_points.iterrows():
    idx.insert(id, well.geometry.bounds)
    
# get the one and only polygon from the district dataset
polygon = gulu_district.geometry.iloc[0]

# get the indexes of wells that intersect bounds of the district
possible_matches_index = list(idx.intersection(polygon.bounds))

# use those indexes to extract the possible matches from the GeoDataFrame
possible_matches = water_points.iloc[possible_matches_index]

# then search the possible matches for precise matches using the slower but more precise method
precise_matches = possible_matches.loc[possible_matches.within(polygon)]

# rebuild the spatial index using the new, smaller dataset
idx = index.Index()
for id, well in precise_matches.iterrows():
	idx.insert(id, well.geometry.bounds)

# declare array to store distances
distances = []

# loop through each water source
for id, house in pop_points.iterrows():

	# use the spatial index to get the index of the closest well
	nearest_well_index = list(idx.nearest(house.geometry.bounds, 1))[0]

	# use the spatial index to get the closest well object from the original dataset
	nearest_well = water_points.iloc[nearest_well_index]

	# store the distance to the nearest well
	distances.append(distance(house.geometry.bounds[0], house.geometry.bounds[1],
		nearest_well.geometry.bounds[0], nearest_well.geometry.bounds[1]))