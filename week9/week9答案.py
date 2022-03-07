"""
Understanding GIS: Practical 8
@author jonnyhuck

Design a footpath between Keswick and the summit of Helvellyn based upon walking efficiency

References:
	https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.contour.html#matplotlib.axes.Axes.contour
	https://rasterio.readthedocs.io/en/latest/api/rasterio.plot.html
	https://rasterio.readthedocs.io/en/stable/topics/plotting.html
	https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html

New Topics:
    Custom Exceptions
"""

from sys import exit
from numpy import linspace
from geopandas import GeoSeries
from math import floor, ceil, hypot
from rasterio import open as rio_open
from shapely.geometry import LineString
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from rasterio.plot import show as rio_show
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import subplots, savefig, get_cmap, colorbar


# this is a very basic exception class that we will use
class PathNotFound(Exception):
    pass


def distance(o, d):
	"""
	* Cartesian distance using Pythagoras
	"""
	return hypot(o[0] - d[0], o[1] - d[1])


def score(cur, nei, data, start, end):
	"""
	* Provides a weighted score by comparing the current node with a neighbouring node.
	* Loosely based on the Nisson score concept: f=g+h
	* In this case, the "h" (heuristic) value is the elevation value of each node.
	"""
	# init score variable
	score = 0

	# current node
	cur_h = data[cur] 				# elevation
	cur_g = distance(cur, end) 		# dist from end
	cur_d = distance(cur, start) 	# dist from start

	# neighbour node
	node_h = data[nei]				# elevation
	node_g = distance(nei, end) 	# dist from end
	node_d = distance(nei, start) 	# dist from start

	# if it is downhill, award the difference in elevation as a score
	if node_h < cur_h:
		score += cur_h - node_h

	# if it is uphill, penalise the score by the difference in elevation
	elif node_h > cur_h:
		score -= node_h - cur_h

	# if it is closer to the end, bonus points
	if node_g < cur_g:
		score += 50

	# if it is further from the start, bonus points
	if node_d > cur_d:
		score += 25

	# return the score
	return score


def a_star(start, end, data):
	"""
	* A-star search algorithm. Moves through nodes in a network (or grid), scores each node's
	* neighbours, and goes to the node with the best score until it finds the end.
	"""

	# initialise sets to handle the data
	closed_cells = set()	# closed set of nodes to avoid
	open_cells = set()		# open set of nodes to evaluate
	path = []				# output list of path nodes

	# add the starting point to begin processing
	open_cells.add(start)

	# keep looping until openCells is empty (we run out of data)
	while open_cells:

		# remove the the next node from openCells
		cur = open_cells.pop()

		# return if we're at the end
		if cur == end:
			return path

		# close off this node to future processing
		closed_cells.add(cur)

		# add the current node to the path
		path.append(cur)

		# make an empty list to hold neighbouring nodes for processing
		options = []

		# load the 8 neighbours to the options list (if they all exist)
		x1, y1 = cur
		if y1 > 0:									# up
			options.append((x1, y1-1))
		if y1 < dem.height-1:						# down
			options.append((x1, y1+1))
		if x1 > 0:									# left
			options.append((x1-1, y1))
		if x1 < dem.width-1:						# right
			options.append((x1+1, y1))
		if x1 > 0 and y1 > 0:						# top left
			options.append((x1-1, y1-1))
		if y1 < dem.height-1 and x1 < dem.width-1:	# bottom right
			options.append((x1+1, y1+1))
		if y1 < dem.height-1 and x1 > 0:			# bottom left
			options.append((x1-1, y1+1))
		if y1 > 0 and x1 <  dem.width-1:			# top right
			options.append((x1+1, y1-1))

		# if the end is a neighbour, add to path and return
		if end in options:
			path.append(end)
			return path

		# see which is the best (biggest) of the scores for the neighbours
		best_score = 0
		for option in options:

			# skip the cell if it has already been assessed
			if option not in closed_cells:

				# get the score and compare it to the best known
				option_score = score(cur, option, data, start, end)

				# if it is better, replace the current best with it
				if option_score > best_score:
					best = option
					best_score = option_score
				else:
					# if not, seal it off so that it is never assessed again
					closed_cells.add(option)

		# add the best node to the open set
		open_cells.add(best)

	# if it gets to here, something went wrong - return an empty list
	raise PathNotFound


# open the specified data file
with rio_open('../data/helvellyn/Helvellyn.tif') as dem:
	# print(dem.profile)

	# read in the data from band 1 (the only band in a DEM)
	band_1 = dem.read(1)

	# set the origin and destination for your path
	origin = dem.index(326794, 523492)		# Keswick town centre
	destination = dem.index(334241, 515107)	# Helvellyn

	# calculate the least cost path using the A* algorithm
	try:
		path = a_star(origin, destination, band_1)

	# catch exception meaning that we have found a path
	except PathNotFound:
		print("Sorry, could not find a path")
		exit()

	# loop through the result and convert it back to coordinates
	coords = [ dem.xy(p[0], p[1]) for p in path ]

	# store as Shapely LineString
	result = LineString(coords)

	# output image
	fig, my_ax = subplots(1, 1, figsize=(16, 10))
	my_ax.set(title="Least Cost Path from Keswick to Helvellyn Summit")

	# set bounds
	buffer = 1000
	my_ax.set_xlim([result.bounds[0] - buffer, result.bounds[2] + buffer])
	my_ax.set_ylim([result.bounds[1] - buffer, result.bounds[3] + buffer])

	# create customised colour map (clipped off the very ends to make it look better)
	terrain = get_cmap('terrain')
	terrain2 = LinearSegmentedColormap.from_list('terrain2', terrain(linspace(0.25, 0.95, 100)))
	
	# draw dem
	rio_show(
		band_1,
		ax=my_ax,
		transform = dem.transform,
		cmap = terrain2,
		)

	# draw dem as contours
	rio_show(
		band_1,
		ax=my_ax,
		contour=True,
		transform = dem.transform,
		colors = ['white'],
		linewidths = [0.5],
		)

	# draw our route as a thick white line then a thinner black line, to give the impression
	#  of a white case around the line to make it stand out better
	GeoSeries(result).plot(
		ax = my_ax,
		linewidth = 5,
		edgecolor = 'white',
		)
	GeoSeries(result).plot(
		ax = my_ax,
		linewidth = 3,
		edgecolor = 'black',
		)

	# add a colour bar
	fig.colorbar(ScalarMappable(norm=Normalize(vmin=floor(band_1.min()), vmax=ceil(band_1.max())), cmap=terrain2), ax=my_ax)

	# add north arrow
	x, y, arrow_length = 0.97, 0.99, 0.1
	my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
		arrowprops=dict(facecolor='black', width=5, headwidth=15),
		ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

	# add scalebar
	my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower right"))

	# save the result
	savefig('./out/8.png', bbox_inches='tight')
	print("done!")
