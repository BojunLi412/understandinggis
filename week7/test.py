# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 10:40:59 2021

@author: 81465
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 21:17:11 2021

@author: Bojun Li
"""
from time import time

# set start time
start_time = time()	# NO CODE ABOVE HERE

from rtree import index
from osmnx import graph_from_file
from networkx import astar_path, NetworkXNoPath
from pyproj import Geod
from sys import exit
from shapely.geometry import LineString
from matplotlib.patches import Patch
from geopandas import read_file, GeoSeries
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import subplots, savefig, Line2D
from pathlib import Path
from pickle import dump, load, HIGHEST_PROTOCOL

def ellipsoidalDistance(a, b):
	"""
	Calculate the 'as the crow flies' distance between two locations along an ellipsoid using
	the Inverse Vincenty method, via the PyProj library.
	"""

	# set which ellipsoid you would like to use
	g = Geod(ellps='WGS84') # this one is a pretty safe bet for global stuff

	# extract nodes
	start = graph.nodes(data=True)[a]
	end = graph.nodes(data=True)[b]

	# compute forward and back azimuths, plus distance
	azF, azB, distance = g.inv(start['x'], start['y'], end['x'], end['y'])

	# we are only interested in the distance
	return distance
'''
# if no pickle data are available
if not Path('../data/osm.pkl').is_file():
	print("No pickle file found, loading data (takes a long time).")

	# load the data into a Python object from the osm file
	graph = graph_from_file('./data/keswick.osm', name="Keswick")

	# pickle the results for future use
	with open('../data/osm.pkl', 'wb') as output:
		dump(graph, output, HIGHEST_PROTOCOL)

else:
	print("Pickle file found successfully, loading data.")

	# extract data from pickle file
	with open('../data/osm.pkl', 'rb') as input:
		graph = load(input)

'''
# create a MultiDiGraph from XML dataset
graph1 = graph_from_file('../data/kaliningrad/keswick.osm', name="Keswick")
'''
# create spatial index from graph
idx = index.Index()
for id, data in list(graph.nodes(data=True)):
	idx.insert(int(id), (data['x'], data['y'], data['x'], data['y']))
# print(idx.count(idx.bounds), len(graph.nodes))

# calculate the 'from' and 'to' node as the nearest to the specified coordinates
fromNode = list(idx.nearest((20.483322,54.692934,20.483322,54.692934), 1))[0]
toNode = list(idx.nearest((20.544863,54.723827,20.544863,54.723827), 1))[0]
# print the 'from' and 'to' nodes to the console
# print(graph.nodes()[fromNode])
# print(graph.nodes()[toNode])
# a=ellipsoidalDistance(fromNode,toNode)
# print(a)
#shortestPath = astar_path(graph, source=fromNode, target=toNode, heuristic=ellipsoidalDistance)
# print(shortestPath)

# use try statement to catch exceptions
try:
	# calculate the shortest path across the network and extract from graph
	shortestPath = astar_path(graph, source=fromNode, target=toNode, heuristic=ellipsoidalDistance)
# 	print(shortestPath)


# catch exception for no path available and exit the script
except NetworkXNoPath:
	print("Sorry, there is no path between those locations in the provided network")
	exit()
#print(shortestPath.length)

# loop through each node in the shortest path and load into list
line = []

for id in shortestPath:
	
	# get the relevant node from the graph with lat lng data
	node = graph.nodes(data=True)[id]

	# load the lat lng data into the lineString
	line.append([node['x'], node['y']])
    
#print(len(line))

# store as a LineString
lineString = LineString(line)

# print(len(lineString.coords))

# convert linestring to GeoSeries and project to UTM zone 34
utm34 = "+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
path = GeoSeries(lineString, crs="+proj=longlat +datum=WGS84 +no_defs").to_crs(utm34)

# open shapefiles for map
buildings = read_file('../data/kaliningrad/buildings.shp').to_crs(utm34)
water = read_file('../data/kaliningrad/water.shp').to_crs(utm34)

# create map axis object, remove axes, set title
fig, my_ax = subplots(1, 1, figsize=(16, 10))
my_ax.axis('off')
my_ax.set(title="The 4 Bridges of Kaliningrad")

# set bounds
buffer = 1000
my_ax.set_xlim([path.geometry.iloc[0].bounds[0] - buffer, path.geometry.iloc[0].bounds[2] + buffer])
my_ax.set_ylim([path.geometry.iloc[0].bounds[1] - buffer, path.geometry.iloc[0].bounds[3] + buffer])

# add the water
water.plot(
	ax=my_ax,
	color='#a6cee3',
	linewidth = 1,
	)

# add the buildings
buildings.plot(
	ax=my_ax,
	color='grey',
	linewidth = 1,
	)

# add the path
path.plot(
	ax=my_ax,
	color='red',
	linewidth = 2,
	)

# manually draw a legend
my_ax.legend([
    Patch(facecolor='grey', label='Buildings'),
    Patch(facecolor='#a6cee3', edgecolor='#a6cee3', label='Water'),
	Line2D([0], [0], color='red', lw=2)],
	['Buildings', 'Water', 'Path'], loc='lower left')

# add north arrow
x, y, arrow_length = 0.99, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=5, headwidth=15),
	ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower right"))

# save the result
savefig(f'out/7.png', bbox_inches='tight')
print("done!")

# report runtime
print(f"completed in: {time() - start_time} seconds")	# NO CODE BELOW HERE      
'''