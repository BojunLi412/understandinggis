# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:12:22 2021

@author: Bojun Li
"""
'''
# here is your list of numbers
numbers = [1,2,3,4,5,6,7,8,9,10,15,30]

# this variable will hold your result, start it at 0
total = 0

# loop through the list
for i in range(12):

	# add each number to the total
	total += numbers[i]

# print the result
print(total)
'''

from geopandas import read_file, GeoSeries
from matplotlib.pyplot import subplots, savefig
from pyproj import Geod
from matplotlib_scalebar.scalebar import ScaleBar
# load the shapefile of countries - this gives a table of 12 columns and 246 rows (one per country)
world = read_file("../data/natural-earth/ne_10m_admin_0_countries.shp")

# print a list of all of the columns in the shapefile
#print(world.columns)


# extract usa and mexico into
# this gives a table only 1 row
usa = world[(world.ISO_A3 == 'USA')] 
mexico = world[(world.ISO_A3 == 'MEX')]

#print(usa.head)
#print(mexico.head)

# extract the geometry column from the (one row) table
# this gives a table with only 1 row AND 1 column
usa_column = usa['geometry']
mexico_column = mexico['geometry']
# extract the raw geometry value that is stored in the first (and only) 1x1 table in usa_column and mexico cllumn
usa_geom = usa_column.iloc[0]
mexico_geom = mexico_column.iloc[0]

# calculate intersection using the raw geometries
border = usa_geom.intersection(mexico_geom)


# create map axis object
my_fig, my_ax = subplots(1, 1, figsize=(16, 10))

# remove axes
my_ax.axis('off')

# plot the border
GeoSeries(border).plot(
  ax = my_ax
	)

# save the image
savefig('./out/first-border.png')


# set which ellipsoid you would like to use
g = Geod(ellps='WGS84')

#print(border)
# initialise a variable to hold the cumulative length
cumulative_length = 0

# loop through each segment in the line
for segment in list(border):

	# THIS LINE NEEDS COMPLETING
	azF, azB, distance = g.inv(segment.coords[0][0], segment.coords[0][1],
		segment.coords[1][0], segment.coords[1][1])

	# add the distance to our cumulative total
	cumulative_length += distance
print(cumulative_length)
# open the graticule dataset
graticule = read_file("../data/natural-earth/ne_110m_graticules_5.shp")
# create map axis object
my_fig, my_ax = subplots(1, 1, figsize=(16, 10))

# remove axes
my_ax.axis('off')

# set title
my_ax.set(title= f"Trump's wall would have been {cumulative_length / 1000:.2f} km long.")

# project border
lambert_conic = '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
border = GeoSeries(border, crs=world.crs).to_crs(lambert_conic)

# set bounds (10000m buffer around the border itself, to give us some context)
buffer = 10000
my_ax.set_xlim([border.geometry.iloc[0].bounds[0] - buffer, border.geometry.iloc[0].bounds[2] + buffer])
my_ax.set_ylim([border.geometry.iloc[0].bounds[1] - buffer, border.geometry.iloc[0].bounds[3] + buffer])

# plot data
usa.to_crs(lambert_conic).plot(
    ax = my_ax,
    color = '#f0e0e0',
    edgecolor = '#660000',
    linewidth = 0.5,
    )
mexico.to_crs(lambert_conic).plot(
    ax = my_ax,
    color = '#e0f0e0',
    edgecolor = '#006600',
    linewidth = 0.5,
    )
border.plot(
    ax = my_ax,
    edgecolor = '#21209e',
    linewidth = 2,
    )
graticule.to_crs(lambert_conic).plot(
    ax=my_ax,
    color='grey',
    linewidth = 1,
    )

# save the result
savefig('out/2.png', bbox_inches='tight')
print("done!")

# add north arrow
x, y, arrow_length = 0.98, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=5, headwidth=15),
	ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.25))


