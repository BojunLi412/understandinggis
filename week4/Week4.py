# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 09:42:35 2021

@author: Bojun Li
"""

from shapely.geometry import LineString
from geopandas import read_file, GeoSeries
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import subplots, savefig, Line2D

def get_effective_area(a, b, c):
    """
    Calculate the area of a triangle made from the points a, b and c
    Simple 'secondary school' approach: area = half * base * height
    """
    return abs( (a[0]-c[0]) * (b[1]-a[1]) - (a[0]-b[0]) * (c[1]-a[1]) ) * 0.5

def visvalingam_whyatt(node_list,n_nodes):
    # loop through each node, excluding the end points
    areas = []
    for i in range(1, len(node_list)-1):
    
      # get the effective area
      area = get_effective_area(node_list[i-1], node_list[i],node_list[i+1])
 
      # append the node and effective area to the list
      areas.append({"point": node_list[i], "area": area})
    # add the end points back in
    areas.insert(0, {"point": node_list[0], "area": 0})
    areas.insert(len(areas), {"point": node_list[len(node_list)-1], "area": 0})  
    # take a copy of the list so that we don't edit the original
    nodes = areas.copy()

     
    # keep going until we run out of nodes
    while len(nodes) > n_nodes:
        min_area = float("inf")
        for i in range(1,len(nodes)-1):
            if min_area > nodes[i]['area']: 
                min_area = nodes[i]['area']
                current_point=i
    # remove the current point from the list
        nodes.pop(current_point)
              
# recalculate effective area to the left of the deleted node
        nodes[current_point-1]['area'] = get_effective_area(nodes[current_point-2]['point'],
    			nodes[current_point-1]['point'], nodes[current_point]['point'])		# left

        # if there is a node to the right of the deleted node, recalculate the effective area
        if current_point < len(nodes)-1:
            nodes[current_point]['area'] = get_effective_area(nodes[current_point-1]['point'],
            				nodes[current_point]['point'], nodes[current_point+1]['point'])	# right
    # extract the nodes and return
    return [ node['point'] for node in nodes ]        

# get the proj string definition for British National Grid (OSGB)
osgb = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs"

# open a dataset of all countries in the world
world = read_file("../data/natural-earth/ne_10m_admin_0_countries.shp")

# extract the UK, project, and extract the geometry (a multipolygon)
uk = world[(world.ISO_A3 == 'GBR')].to_crs(osgb).geometry.iloc[0]

# report geometry type
print(f"geometry type: {uk.geom_type}")
coord_list = []
biggestarea = 0
# extract the biggest area with loop
for poly in uk:
    if poly.area > biggestarea: 
        biggestarea = poly.area
        # extract the list of coordinates from the polygon
        coord_list =list(poly.exterior.coords)
print(biggestarea)
# set the percentage of nodes that you want to remove
SIMPLIFICATION_PERC = 98

# how many nodes do we need?
n_nodes = int(len(coord_list) / 100.0 * (100 - SIMPLIFICATION_PERC))

# remove one node and overwrite it with the new, shorter list
simplified_nodes = visvalingam_whyatt(coord_list, n_nodes)

# make a linestring out of the coordinates
before_line = LineString(coord_list)
print(f"original node count: {len(coord_list)}")
print(f"original length: {before_line.length / 1000:.2f}km")
print() # print a blank line

# make the resulting list of coordinates into a line
after_line = LineString(simplified_nodes)
print(f"simplified node count: {len(simplified_nodes)}")
print(f"simplified length: {after_line.length / 1000:.2f}km")
print() # print a blank line

# create map axis object
fig, my_ax = subplots(1, 1, figsize=(16, 10))

# remove axes
my_ax.axis('off')

# set title
my_ax.set(title=f"The Coastline of Mainland Great Britain \n({before_line.length / 1000:.0f}km or {after_line.length / 1000:.0f}km long)")

# add the new coastline
GeoSeries(after_line).plot(
    ax=my_ax,
    color='red',
		alpha=0.5,
    linewidth = 0.6,
	)

# add the original coastline
GeoSeries(before_line).plot(
    ax=my_ax,
    color='blue',
		alpha=0.5,
    linewidth = 0.6,
	)

# manually draw a legend
my_ax.legend([
	Line2D([0], [0], color='blue', lw=0.6),
	Line2D([0], [0], color='red',  lw=0.6)],
	['Original', f'{SIMPLIFICATION_PERC}% Simplified'], loc='lower right')

# add north arrow
x, y, arrow_length = 0.95, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=5, headwidth=15),
	ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.25))

# save the result
savefig(f'out/4.png', bbox_inches='tight')
print("done!")

