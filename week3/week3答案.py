# -*- coding: utf-8 -*-



"""
Understanding GIS: Practical 3
@author jonnyhuck

Use spatial index to calculate the mean distance to a water supply in Gulu District

References:
	https://epsg.io/21096
	http://geopandas.org/reference.html
	https://toblerity.org/rtree/examples.html
	https://toblerity.org/rtree/tutorial.html#query-the-index
	https://geoffboeing.com/2016/10/r-tree-spatial-index-python
	https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html
  https://matplotlib.org/3.3.3/tutorials/intermediate/legend_guide.html

New Topics:
  Creating Functions
	RTree
	Filtering GeoDataFrames with loc & iloc
"""

from sys import exit
from math import sqrt
from rtree import index
from geopandas import read_file
from matplotlib.pyplot import subplots, savefig
from matplotlib_scalebar.scalebar import ScaleBar

def distance(x1, y1, x2, y2):
	"""
	* Use Pythagoras' theorem to measure the distance
	* This is acceptable in this application because the distances are so short
	*  and comutational efficiency is a key consideration
	"""
	return sqrt((x1 - x2)**2 + (y1 - y2)**2)


# read in shapefiles, ensure that they all have the same CRS
pop_points = read_file("../data/gulu/pop_points.shp")
water_points = read_file("../data/gulu/water_points.shp").to_crs(pop_points.crs)
gulu_district = read_file("../data/gulu/district.shp").to_crs(pop_points.crs)


''' STAGE 1 - FILTERING THE DATA WITH A SPATIAL QUERY - EXTRACT ONLY WELLS INSIDE THE DISTRICT '''

# initialise an rtree index
idx = index.Index()

# loop through each row in the wells dataset and construct spatial index
for id, well in water_points.iterrows():
	idx.insert(id, well.geometry.bounds)

# get the one and only polygon from the district dataset (buffer by 10.5km to account for edge cases)
polygon = gulu_district.geometry.iloc[0].buffer(10500)

# report how many wells there are at this stage
print(f"Initial wells: {len(water_points.index)}")

# get the indexes of wells that intersect bounds of the district, then extract extract the possible matches
possible_matches = water_points.iloc[list(idx.intersection(polygon.bounds))]
print(f"Potential wells: {len(possible_matches.index)}")

# then search the possible matches for precise matches using the slower but more precise method
precise_matches = possible_matches.loc[possible_matches.within(polygon)]
print(f"Filtered wells: {len(precise_matches.index)}")


''' STAGE 2 - NEAREST NEIGHBOUR - find distance to nearest well for everybody'''

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

# calculate the mean
mean = round(sum(distances) / len(distances))

# store distance to nearest well
pop_points['nearest_well'] = distances

# output the result
print(f"Minimum distance to water in Gulu District is {round(min(distances))}m.")
print(f"Mean distance to water in Gulu District is {mean}m.")
print(f"Maximum distance to water in Gulu District is {round(max(distances))}m.")

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