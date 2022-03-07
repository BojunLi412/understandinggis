# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:41:32 2021

@author: Bojun Li
"""
# THIS IMPORTS A SQUARE ROOT FUNCTION, WHICH IS A BIG HINT!!

from math import sqrt
from geopandas import read_file
from rtree import index
from matplotlib.pyplot import subplots, savefig
from matplotlib_scalebar.scalebar import ScaleBar

def distance(x1, y1, x2, y2):   
    """
    * Use Pythagoras' theorem to measure the distance
    * This is acceptable in this application because the distances are so short
    *  and computational efficiency is a key consideration
    """
    return sqrt((x1-x2)**2+(y1-y2)**2)


'''
a= distance(1, 1, 4, 5)
print(a)
'''

# read in shapefiles, ensure that they all have the same CRS
pop_points = read_file("../data/gulu/pop_points.shp")
water_points = read_file("../data/gulu/water_points.shp").to_crs(pop_points.crs)
gulu_district = read_file("../data/gulu/district.shp").to_crs(pop_points.crs)



'''
print(pop_points.crs.to_epsg())
print()
print(water_points.crs.to_epsg())
print()
print(gulu_district.crs.to_epsg())

print(f"population points: {len(pop_points.index)}")
print(f"Initial wells: {len(water_points.index)}")
'''
# initialise an rtree index
idx = index.Index()
print(water_points.iterrows)
'''
# loop through each row and construct spatial index on the larger dataset (population)
for id, well in water_points.iterrows():
    idx.insert(id, well.geometry.bounds)
# get the one and only polygon from the district dataset
polygon = gulu_district.geometry.iloc[0].buffer(10500)

# how many rows are we starting with?
print(f"Initial wells: {len(water_points.index)}")

# get the indexes of wells that intersect bounds of the district
possible_matches_index = list(idx.intersection(polygon.bounds))

# use those indexes to extract the possible matches from the GeoDataFrame
possible_matches = water_points.iloc[possible_matches_index]

# how many rows are left now? 
print(f"Filtered wells: {len(possible_matches.index)}")

# then search the possible matches for precise matches using the slower but more precise method
precise_matches = possible_matches.loc[possible_matches.within(polygon)]



# how many rows are left now?
print(f"Filtered wells: {len(precise_matches.index)}")
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
mean = round(sum(distances) / len(distances))
print(len(distances))
print(distances[:5])
# store distance to nearest well
pop_points['nearest_well'] = distances
print(f"Maximum distance to water in Gulu District is {round(max(distances))}m.")
print(f"Mean distance to water in Gulu District is {mean}m.")
print(f"Minimum distance to water in Gulu District is {round(min(distances))}m.")
# create map axis object
fig, my_ax = subplots(1, 1, figsize=(16, 10))

# remove axes
my_ax.axis('off')

# add title
my_ax.set(title="Distance to Nearest Well, Gulu District, Uganda")

# add the district boundary
gulu_district.plot(
    ax = my_ax,
    color = 'white',
    linewidth = 1,
	edgecolor = 'black',
    )

# plot the locations, coloured by distance to water
pop_points.plot(
    ax = my_ax,
    column = 'nearest_well',
    linewidth = 0,
	markersize = 1,
    cmap = 'RdYlBu_r',
    scheme = 'quantiles',
    legend = 'True',
    legend_kwds = {
        'loc': 'lower right',
        'title': 'Distance to Nearest Well'
        }
    )

# add north arrow
x, y, arrow_length = 0.98, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=5, headwidth=15),
	ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.25))

# save the result
savefig('out/3.png', bbox_inches='tight')
print("done!")
'''
