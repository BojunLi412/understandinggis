# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:57:12 2021

@author: Bojun Li
"""


from geopandas import read_file
from pyproj import Geod, CRS, Transformer
from numpy.random import uniform
from numpy import arange
from shapely.geometry import Polygon
from math import sin, cos, radians, hypot


def offset(x, y, distance, direction):
    """
    * Offset a location by a given distance and direction
    """
    x2 = x + cos(radians(direction)) * distance
    y2 = y + sin(radians(direction)) * distance
    return (x2, y2)


def evaluate_distortion(transformer, minx, miny, maxx, maxy, sample_number):

    # calculate the required number of random locations (x and y separately) plus radius
    xs = uniform(low=minx, high=maxx, size=sample_number)
    ys = uniform(low=miny, high=maxy, size=sample_number)
    rs = uniform(low=1000, high=1000000, size=sample_number)
    
    # offset distances
    forward_azimuths = arange(0, 360, 22.5)
    n = len(forward_azimuths)
    
    # loop through the points
    planar_areas = []
    shape_indices = []
    ellipsoidal_areas = []
    
    for x, y, r in zip(xs, ys, rs):
        
        # construct a circle around the centre point on the ellipsoid
        lons, lats, bAz = g.fwd([x]*n, [y]*n, forward_azimuths, [r]*n)
    
        # project the result, calculate area, append to the list
        e_coords = [ transformer.transform(lon, lat, direction='FORWARD') for lon, lat in zip(lons, lats)]
        ellipsoidal_areas.append(Polygon(e_coords).area) 
        # transform the centre point to the projected CRS
        px, py = transformer.transform(x, y, direction='FORWARD')
        
        # construct a circle around the projected point on a plane, calculate area, append to list
        p_coords = [ offset(px, py, r, az) for az in forward_azimuths ]
        planar_areas.append(Polygon(p_coords).area)
        # get radial distances from the centre to each of the 16 points on the circle
        ellipsoidal_radial_distances = [ hypot(px - ex, py - ey) for ex, ey in e_coords ]
        
        # get the sum of the distances, and the expected value for each distance
        total_radial_dist = sum(ellipsoidal_radial_distances)
        expected_distance = total_radial_dist / n
        # get the difference between the actual and expected radial distance for each 'spoke'
        shape_distortion = [ abs((expected_distance / total_radial_dist) - (d / total_radial_dist)) for d in ellipsoidal_radial_distances ]
        shape_indices.append(sum(shape_distortion))
    # calculate shape distortion
    Es = sum(shape_indices) / len(shape_indices)    
    # get the sum of the differences between the ellipsoidal and planar circle areas
    diff_sum = 0
    for e, p in zip(ellipsoidal_areas, planar_areas):
      diff_sum += abs(e - p) / abs(e + p)
    
    # calculate the Ea and Ka numbers
    Ea = 1 / sample_number * diff_sum
    Ka = (1 + Ea) / (1 - Ea)
    return Ep, Es, Ea, Ka


world = read_file("../data/natural-earth/ne_10m_admin_0_countries.shp")

# extract usa into
# this gives a table only 1 row

country = world[(world.ISO_A3 == 'GRL')] 

# get the bounds of the country
minx, miny, maxx, maxy = country.total_bounds

# set which ellipsoid you would like to use
g = Geod(ellps='WGS84')

proj_string  = '+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs'

geo_string = "+proj=longlat +datum=WGS84 +no_defs "

# initialise a PyProj Transformer to transform coordinates
transformer = Transformer.from_crs(CRS.from_proj4(geo_string), CRS.from_proj4(proj_string), always_xy=True)

'''
Function
'''

sample_number=1000
# calculate the required number of random locations (x and y separately) plus radius
xs = uniform(low=minx, high=maxx, size=sample_number)
ys = uniform(low=miny, high=maxy, size=sample_number)
rs = uniform(low=1000, high=1000000, size=sample_number)

# offset distances
forward_azimuths = arange(0, 360, 22.5)
n = len(forward_azimuths)

# loop through the points
planar_areas = []
shape_indices = []
ellipsoidal_areas = []

for x, y, r in zip(xs, ys, rs):
    
    # construct a circle around the centre point on the ellipsoid
    lons, lats, bAz = g.fwd([x]*n, [y]*n, forward_azimuths, [r]*n)

    # project the result, calculate area, append to the list
    e_coords = [ transformer.transform(lon, lat, direction='FORWARD') for lon, lat in zip(lons, lats)]
    ellipsoidal_areas.append(Polygon(e_coords).area) 
    # transform the centre point to the projected CRS
    px, py = transformer.transform(x, y, direction='FORWARD')
    
    # construct a circle around the projected point on a plane, calculate area, append to list
    p_coords = [ offset(px, py, r, az) for az in forward_azimuths ]
    planar_areas.append(Polygon(p_coords).area)
    # get radial distances from the centre to each of the 16 points on the circle
    ellipsoidal_radial_distances = [ hypot(px - ex, py - ey) for ex, ey in e_coords ]
    
    # get the sum of the distances, and the expected value for each distance
    total_radial_dist = sum(ellipsoidal_radial_distances)
    expected_distance = total_radial_dist / n
    # get the difference between the actual and expected radial distance for each 'spoke'
    shape_distortion = [ abs((expected_distance / total_radial_dist) - (d / total_radial_dist)) for d in ellipsoidal_radial_distances ]
    shape_indices.append(sum(shape_distortion))
# calculate shape distortion
Es = sum(shape_indices) / len(shape_indices)    
# get the sum of the differences between the ellipsoidal and planar circle areas
diff_sum = 0
for e, p in zip(ellipsoidal_areas, planar_areas):
  diff_sum += abs(e - p) / abs(e + p)

# calculate the Ea and Ka numbers
Ea = 1 / sample_number * diff_sum
Ka = (1 + Ea) / (1 - Ea)
