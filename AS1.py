"""
Understanding GIS: Assessment 1
@author 10621119

Calculate the length of the World's Shortest Border, as fast as possible
"""

from time import time

# set start time
start_time = time()	# NO CODE ABOVE HERE

'''ALL CODE MUST BE INSIDE HERE'''

'''FIND THE MIN BORDER AND INTERSECTION COUNTRIES'''

from geopandas import read_file,GeoSeries
from rtree import index
from pyproj import Geod
from matplotlib.pyplot import subplots, savefig
from matplotlib_scalebar.scalebar import ScaleBar

# set which ellipsoid you would like to use
g = Geod(ellps='WGS84')

# load the shapefile of countries and the graticule
world = read_file("./data/ne_10m_admin_0_countries.shp")
graticule = read_file("./data/ne_110m_graticules_5.shp")

# load all countries into spatisal index
idx = index.Index()

# Initialise mini_borderlength variables as the max length 
mini_borderlength = float('inf')

# Initialise list for final intersection countries
final_countriess = []

# Initialise the minimum of border
min_border = []

# Initialise the checked index 
checked = []

# loop through each row of data, save the checked data in a list and conver  country list to geometry
for id1, country_A in world.iterrows():
    
    # load all countries into spatisal index 
    idx.insert(id1, country_A.geometry.bounds)
    
    # save the index number that has been run
    checked.append(index)
    
    # convert list countries to geometry
    geoms = [country_A.geometry] 

    # loop the list of geometries and get intersections
    for country_B in geoms:
        
        # get the index of the intersected countries and exclude the former checked data to avoid redundant calculation
        possible_matches_index = [ x for x in list(idx.intersection(country_B.bounds)) if x not in checked ]
        
        # skip not match countries
        if len(possible_matches_index) != 0:
          
            # get a GeoDataFrame for potential matches
            possible_matches = world.iloc[possible_matches_index]
    
            # use the spatial index to get the match items from the possible_matches,it is slower but more precise method
            precise_matches = possible_matches.loc[possible_matches.intersects(country_B)]
    
            # loop through precise matches
            for id2, country_C in precise_matches.iterrows():
    
                # calculate the intersection
                border = country_A.geometry.intersection(country_C.geometry)
               
                # if it returns as any other types it can be ignored as is short
                if(border.type == "MultiLineString"):
    
                    # initialise a variable to hold the cumulative length
                    cumulativeLength = 0
    
                    # loop through each segment in the line
                    for segment in list(border):
    
                        # calculate the forwrda azimuth, backward azimuth and direction of the current segment
                        azF, azB, distance = g.inv(segment.coords[0][0], segment.coords[0][1],
                            segment.coords[1][0], segment.coords[1][1])
    
                        # add the distance to our cumulative total
                        cumulativeLength += distance
                      
                    # find the shortest border length
                    if cumulativeLength < mini_borderlength:
                        
                        # Store the shortest border 
                        min_border = border
                        
                        # Store the shortest border length
                        mini_borderlength = cumulativeLength
                        
                        # Store the two countries with the shortest border
                        final_countries = [country_A, country_C]
         
print(f"The shortest border is between {final_countries[0].NAME} and {final_countries[1].NAME}, at {mini_borderlength:.2f}metres.")

#Store the ISO_A3 values of the two countries with the shortest border
Reuslt_Country_1= world.loc[(world.ISO_A3 == final_countries [0].ISO_A3)]
Reuslt_Country_2= world.loc[(world.ISO_A3 == final_countries [1].ISO_A3)]    

''' EXPORT THE RESULT TO THE MAP '''

# create map axis object
my_fig, my_ax = subplots(1, 1, figsize=(16, 10))

# remove axes
my_ax.axis('off')

# set title
my_ax.set(title = f"The shortest border is between {final_countries[0].NAME} and {final_countries[1].NAME}, at {mini_borderlength:.2f}metres.")


# change the min border project 
min_border = GeoSeries(min_border,crs=world.crs).to_crs('+proj=eqearth')

# set bounds (50m buffer around the border itself, to give us some context)
buffer = 50
my_ax.set_xlim([min_border.geometry.iloc[0].bounds[0] - buffer, min_border.geometry.iloc[0].bounds[2] + buffer])
my_ax.set_ylim([min_border.geometry.iloc[0].bounds[1] - buffer, min_border.geometry.iloc[0].bounds[3] + buffer])

# plot data
Reuslt_Country_1.to_crs('+proj=eqearth').plot(
    ax = my_ax,
    color = '#f0e0e0',
    edgecolor = '#660000',
    linewidth = 0.5,
    )
Reuslt_Country_2.to_crs('+proj=eqearth').plot(
    ax = my_ax,
    color = '#e0f0e0',
    edgecolor = '#006600',
    linewidth = 0.5,
    )
min_border.plot(
    ax = my_ax,
    edgecolor = '#21209e',
    linewidth = 2,
    )
graticule.to_crs('+proj=eqearth').plot(
    ax=my_ax,
    color='grey',
    linewidth = 1,
    )

# add north arrow
x, y, arrow_length = 0.98, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=5, headwidth=15),
	ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.25))

# save the result
savefig('out/shortest_border.png', bbox_inches='tight')

# report runtime
print(f"completed in: {time() - start_time} seconds")	# NO CODE BELOW HERE      
