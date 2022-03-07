from geopandas import read_file

# load the countries data into a GeoDataFrame
world = read_file("../data/natural-earth/ne_50m_admin_0_countries.shp")
graticule = read_file("../data/natural-earth/ne_110m_graticules_15.shp")
bbox = read_file("../data/natural-earth/ne_110m_wgs84_bounding_box.shp")

# reproject all three layers to equal earth
world = world.to_crs('+proj=eqearth')
graticule = graticule.to_crs('+proj=eqearth')
bbox = bbox.to_crs('+proj=eqearth')

''' Output Population Density '''

world['pop_density'] = world.POP_EST / (world.area / 1000000)
# test the contents
# print(world.head())
from matplotlib.pyplot import subplots, savefig
# create map axis object
my_fig, my_ax = subplots(1, 1, figsize=(16, 10))

# add bounding box and graticule layers
bbox.plot(
    ax = my_ax,
    color = 'lightgray',
    linewidth = 0,
    )

# plot the countries
world.plot(
    ax = my_ax,
    column = 'pop_density',
    linewidth = 0.5,
    edgecolor = 'gray',
    cmap = 'OrRd',
    scheme = 'quantiles',
    legend = 'True',
    legend_kwds = {
        'loc': 'lower left',
        'title': 'Population Density'
    }
    )

# plot the graticule
graticule.plot(
    ax = my_ax,
    color = 'white',
    linewidth = 0.5,
    )

# turn off the visible axes on the map
my_ax.axis('off')

# add title 
my_ax.set(title="Population Density: Equal Earth Coordinate Reference System")


# save the result
savefig('./out/1.png')
print("done!")

#print(world.columns)

''' Output GDP Per Capita '''
my_fig, my_ax1= subplots(1, 1, figsize=(16, 10))
world['GDP_Per_Capita'] = world.GDP_MD_EST / world.POP_EST

bbox.plot(
    ax = my_ax1,
    color = 'lightgray',
    linewidth = 0,
    )

# plot the countries
world.plot(
    ax = my_ax1,
    column = 'GDP_Per_Capita',
    linewidth = 0.5,
    edgecolor = 'gray',
    cmap = 'OrRd',
    scheme = 'quantiles',
    legend = 'True',
    legend_kwds = {
        'loc': 'lower left',
        'title': 'GDP Per Capita'
    }
    )

# plot the graticule
graticule.plot(
    ax = my_ax1,
    color = 'white',
    linewidth = 0.5,
    )

# turn off the visible axes on the map
my_ax1.axis('off')

# add title 
my_ax1.set(title="GDP Per Capita: Equal Earth Coordinate Reference System")


# save the result
savefig('./out/2.png', bbox_inches='tight')
