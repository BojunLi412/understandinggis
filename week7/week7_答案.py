"""
Understanding GIS: Practical 7
@author jonnyhuck

Find the shortest path between two points in Kaliningrad

References:
	https://geoffboeing.com/2016/11/osmnx-python-street-networks/
	https://networkx.github.io/documentation/stable/reference/algorithms/euler.html
	https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.shortest_paths.astar.astar_path.html#networkx.algorithms.shortest_paths.astar.astar_path
	https://github.com/gboeing/osmnx/issues/247
	https://toblerity.org/rtree/

New Topics:
    Exceptions
"""
from sys import exit
from pyproj import Geod
from rtree import index
from pathlib import Path
from osmnx import graph_from_file
from matplotlib.patches import Patch
from shapely.geometry import LineString
from geopandas import read_file, GeoSeries
from networkx import astar_path, NetworkXNoPath
from pickle import dump, load, HIGHEST_PROTOCOL
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import subplots, savefig, Line2D

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


# if no pickle data are available
if not Path('../data/osm.pkl').is_file():
	print("No pickle file found, loading data (takes a long time).")

	# load the data into a MultiDiGraph from the osm file
	graph = graph_from_file('../data/kaliningrad/kaliningrad.xml', name="Kaliningrad")

	# pickle the results for future use
	with open('../data/osm.pkl', 'wb') as output:
		dump(graph, output, HIGHEST_PROTOCOL)
      
else:
	print("Pickle file found successfully, loading data.")

	# extract data from pickle file
	with open('../data/osm.pkl', 'rb') as input:
		graph = load(input)

# create spatial index from graph
idx = index.Index()
for id, data in list(graph.nodes(data=True)):
	idx.insert(int(id), (data['x'], data['y'], data['x'], data['y']))
# print(idx.count(idx.bounds), len(graph.nodes))

# calculate the 'from' and 'to' node as the nearest to the specified coordinates
fromNode = list(idx.nearest((20.483322,54.692934,20.483322,54.692934), 1))[0]
toNode = list(idx.nearest((20.544863,54.723827,20.544863,54.723827), 1))[0]
# print(graph.nodes()[fromNode], graph.nodes()[toNode])
# use try statement to catch exceptions
try:
	# calculate the shortest path across the network and extract from graph
	shortestPath = astar_path(graph, source=fromNode, target=toNode, heuristic=ellipsoidalDistance)

# catch exception for no path available
except NetworkXNoPath:
	print("Sorry, there is no path between those locations in the provided network")
	exit()

# loop through each node in the shortest path and load into list
line = []
for n in shortestPath:

	# get the relevant node from the graph with lat lng data
	node = graph.nodes(data=True)[n]

	# load the lat lng data into the lineString
	line.append([node['x'], node['y']])

# store as a LineString
lineString = LineString(line)

# convert linestring to GeoSeries and project to UTM zone 34
utm34 = "+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
path = GeoSeries(lineString, crs="+proj=longlat +datum=WGS84 +no_defs").to_crs(utm34)

# open buildings data (use pickle if available)
if not Path('../data/buildings.pkl').is_file():
	buildings = read_file('../data/kaliningrad/buildings.shp').to_crs(utm34)
	with open('../data/buildings.pkl', 'wb') as output:
		dump(buildings, output, HIGHEST_PROTOCOL)
else:
	with open('../data/buildings.pkl', 'rb') as input:
		buildings = load(input)

# open water data (use pickle if available)
if not Path('../data/water.pkl').is_file():
	water = read_file('../data/kaliningrad/water.shp').to_crs(utm34)
	with open('../data/water.pkl', 'wb') as output:
		dump(water, output, HIGHEST_PROTOCOL)
else:
	with open('../data/water.pkl', 'rb') as input:
		water = load(input)

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
savefig(f'./out/7.png', bbox_inches='tight')
print("done!")