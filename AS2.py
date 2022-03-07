"""
Understanding GIS: Assessment 2
@author 10621119

An Implementation Weighted Redistribution Algorithm (Huck et al.)
"""
from time import time

# set start time
start_time = time()    # NO CODE ABOVE HERE

'''ALL CODE MUST BE INSIDE HERE'''

from rasterio import open as rio_open
from geopandas import read_file,GeoSeries
from rtree import index
from math import sqrt,pi,floor,ceil
from shapely.geometry import Point
from numpy import zeros,column_stack
from numpy.random import uniform
from rasterio.plot import show as rio_show
from matplotlib.pyplot import subplots, savefig,get_cmap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from skimage.draw import disk

# crate a function to get the largest weighting value point in random points
def get_largest_weighting_value_point(randompoints):
        
    # defination the max_value as the 0
    max_value = 0
    
    # loop every random points's raster value to gain the max value point
    for randompoint in randompoints:
        
        # convert coordinates of random points to image space
        row,col = weightingSurface.index(randompoint.x, randompoint.y)  
        
        # get the value of the raster data
        value = band[row][col]
        
        # judge whether it is the point with the largest value, if it is, save the largest value and point
        if value > max_value:
            max_value = value
            max_point = randompoint
            
    # return the max_value point
    return(max_point)      

# crate a function to generate a certain number of random points within the 
def generate_random_points(admin,sample_number):
    
    # get boundary of polygon's maximum value of the endpoint coordinates in the four directions
    min_x,min_y,max_x,max_y=admin.geometry.bounds
    
    # creat an empty list to store the within admin's random points
    random_points_within=[]
    
    # loop each tweets in the adminarea
    for i in range(sample_number):
        
        random_point = Point([uniform(min_x,max_x) , uniform(min_y,max_y)])
        
        #Check if the point is within the admin
        if random_point.within(admin.geometry):
            random_points_within.append(random_point)
            
    return random_points_within


# get a function to calculate the distance between two points
def distance(x1,x2,y1,y2):
    return sqrt((x1 - x2)**2 + (y1 - y2)**2)

# get the proj string definition for British National Grid (OSGB)
osgb = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs"

#load the level3-tweets point data
tweets = read_file("./data/wr/level3-tweets-subset.shp").to_crs(osgb)

#load the great manchester polygon data 
admins = read_file("./data/wr/gm-districts.shp").to_crs(osgb)

#[user defined] desired influence of weighting surface
w = 20

#[user defined] desired level of spatial ambiguity
s = 0.1

# initialise an rtree index
idx = index.Index()

#Create a spatial index of tweets
for id, tweet in tweets.iterrows(): 
    idx.insert(id, tweet.geometry.bounds)

# open the population raster data
with rio_open('./data/wr/100m_pop_2019.tif') as weightingSurface:

    # read in the pop data  (the only band in this raster)
    band = weightingSurface.read(1)

    # create a new 'band' of raster data the same size
    output = zeros((weightingSurface.height, weightingSurface.width))

    # loop through each admin in admins
    for id,admin in admins.iterrows():

        # calculate r based on Equation 1
        r = sqrt((admin.geometry.area * s) / pi)

        #convert r to pixels
        r_pixels = int(r/weightingSurface.res[0])

        try:

            # find out the existing tweets points in each administrative district
            possible_tweets = tweets.iloc[list(idx.intersection(admin.geometry.bounds))]
            precise_tweets = possible_tweets.loc[possible_tweets.within(admin.geometry)] 

            # loop every tweet which is in the admin
            for i in range(len(precise_tweets)):

                # generate random points within the admin
                randompoints = generate_random_points(admin,w)

                # get the largest weighting value of point as the seed point
                seed = get_largest_weighting_value_point(randompoints)

                # convert coordinates of seed to image space
                row_seed,col_seed=weightingSurface.index(seed.x,seed.y)
                
                # for row, col in row_seed,col_seed:
                for row, col in column_stack(disk((row_seed, col_seed), r_pixels)):
                    
                    # use equation 2 to get the output value
                    output[row][col] += 1-(sqrt((row_seed-row)**2+(col_seed-col)**2)/r_pixels)

        except Exception:
            pass
        

# create map axis object
fig, my_ax = subplots(1, 1, figsize=(16, 10))
    
# set the title
my_ax.set(title="Twitters of the Royal Wedding in Greater Manchester")

# get colour map
colour_map = get_cmap('YlOrRd')
    
# draw output 
rio_show(
    output,
    ax = my_ax,
    transform = weightingSurface.transform,
    cmap = 'YlOrRd'
    )

# draw administrative areas
GeoSeries(admins.geometry.boundary).plot(
    ax = my_ax,
    linewidth = 0.2,
    edgecolor = 'black'
    )

# add north arrow
x, y, arrow_length = 0.97, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y - arrow_length),
                   arrowprops=dict(facecolor='black', width=5, headwidth=15),
                   ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add colorbar
fig.colorbar(ScalarMappable(norm=Normalize(vmin=floor(output.min()), vmax=ceil(output.max())), cmap=colour_map), ax=my_ax)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower right"))

# save the result
savefig('./out/manchester_tweets.png', bbox_inches='tight')

# report runtime
print(f"completed in: {time() - start_time} seconds")    # NO CODE BELOW HERE