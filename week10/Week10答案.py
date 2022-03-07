"""
Understanding GIS: Practical 10
@author jonnyhuck

Calculate a Viewshed using Bresenham's Line and Midpoint Circle Algorithm

References:
	https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.contour.html#matplotlib.axes.Axes.contour
	https://rasterio.readthedocs.io/en/latest/api/rasterio.plot.html
	https://rasterio.readthedocs.io/en/stable/topics/plotting.html
	https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
	http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=How%20Visibility%20works

Library:
	conda install -c conda-forge scikit-image
"""

from sys import exit
from numpy import linspace
from geopandas import GeoSeries
from shapely.geometry import Point
from math import hypot, floor, ceil
from numpy import zeros, column_stack
from rasterio import open as rio_open
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from rasterio.plot import show as rio_show
from skimage.draw import line, circle_perimeter
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import subplots, savefig, get_cmap


def adjust_height(height, distance, earth_diameter=12740000, refraction_coefficient=0.13):
	"""
	* Adjust the apparant height of an object at a certain distance, accounting for the
	* 	curvature of the earth and atmospheric refraction
	"""
	# this is the arcGIS version
	# a = distance**2 / earth_diameter
	# return height - a + refraction_coefficient * a

	# this is the QGIS version
	return height - (distance**2 / earth_diameter) * (1 - refraction_coefficient)


def line_of_sight(r0, c0, height0, r1, c1, height1, radius, dem_data, output):
	"""
	 * Use Bresenham's Line algorithm to calculate a line of sight from one point to another point, 
	 *	returning a list of visible cells
	"""
	# init variable for biggest dydx so far (starts at -infinity)
	max_dydx = -float('inf')

	# loop along the pixels in the line (exclusing the first one)
	for r, c in column_stack(line(r0, c0, r1, c1))[1:]:

		# calculate distance travelled so far (in pixels)
		dx = hypot(c0 - c, r0 - r)
		
		# stop if we have gone too far
		if dx > radius:
			break

		try:
			# calculate the current value for dy / dx
			base_dydx = (adjust_height(dem_data[(r, c)], dx*50) - height0) / dx
			tip_dydx = (adjust_height(dem_data[(r, c)] + height1, dx*50) - height0) / dx
		
		# if we go off the edge of the data, stop looping (EAFP)
		except IndexError:
			break
		
		# if the tip dydx is bigger than the previous max, it is visible
		if (tip_dydx > max_dydx):
			output[(r, c)] = 1

		# if the base dydx is bigger than the previous max, update
		max_dydx = max(max_dydx, base_dydx)

	# return updated output surface
	return output


def viewshed(x0, y0, radius_m, observer_height, target_height, dem_data, d):
	"""
	* Use Midpoint Circle algorithm to determine endpoints for viewshed
	"""

	# convert origin coordinates to image space
	r0, c0 = d.index(x0, y0)

	# convert the radius (m) to pixels
	radius_px = int(radius_m / d.res[0])

	try:
		# get the observer height (above sea level)
		height0 = dem_data[(r0, c0)] + observer_height
	
	# if it is not within the dataset, in form the user and exit (EAFP)
	except IndexError:
		print(f"Sorry, {x0, y0} is not within the elevation dataset.")
		exit()

	# create output array at the same dimensions as data for viewshed
	output = zeros(dem_data.shape)

	# set the start location as visible automatically
	output[(r0, c0)] = 1

	# get pixels in the perimeter of the viewshed
	for r, c in column_stack(circle_perimeter(r0, c0, radius_px*3)):

		# calculate line of sight to each pixel, pass output and get a new one back each time
		output = line_of_sight(r0, c0, height0, r, c, target_height, radius_px, dem_data, output)

	# return the resulting viewshed
	return output


# open the elevation data file
with rio_open("../data/helvellyn/Helvellyn-50.asc") as dtm_input:

	# read the data out of band 1 in the dataset
	dem_data = dtm_input.read(1)

	# set origin for viewshed
	x0, y0 = 332000, 514000

	# calculate the viewshed
	output = viewshed(x0, y0, 5000, 1.8, 100, dem_data, dtm_input)

	# output image
	fig, my_ax = subplots(1, 1, figsize=(16, 10))
	my_ax.set(title="Viewshed Analysis")

	# create customised colour map (clipped off the very ends to make it look better)
	terrain = get_cmap('terrain')
	terrain2 = LinearSegmentedColormap.from_list('terrain2', terrain(linspace(0.25, 0.95, 100)))

	# draw dem
	rio_show(
		dem_data,
		ax=my_ax,
		transform = dtm_input.transform,
		cmap = terrain2,
		)

	# draw dem as contours
	rio_show(
		dem_data,
		ax=my_ax,
		contour=True,
		transform = dtm_input.transform,
		colors = ['white'],
		linewidths = [0.5],
		)

	# add viewshed
	rio_show(
		output,
		ax=my_ax,
		transform=dtm_input.transform,
		cmap = LinearSegmentedColormap.from_list('binary_viewshed', [(0, 0, 0, 0), (1, 0, 0, 0.5)], N=2)
		)

	# add origin point
	GeoSeries(Point(x0, y0)).plot(
		ax = my_ax,
		color = 'black',
		)

	# add a colour bar
	fig.colorbar(ScalarMappable(norm=Normalize(vmin=floor(dem_data.min()), vmax=ceil(dem_data.max())), cmap=terrain2), ax=my_ax)

	# add north arrow
	x, y, arrow_length = 0.97, 0.99, 0.1
	my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
		arrowprops=dict(facecolor='black', width=5, headwidth=15),
		ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

	# add scalebar
	my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower right"))

	# save the result
	savefig('./out/10.png', bbox_inches='tight')
	print("done!")
