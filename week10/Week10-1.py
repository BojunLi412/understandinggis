"""
Created on Mon Apr 26 09:59:21 2021

@author: Bojun Li
"""

from rasterio import open as rio_open
from rasterio.plot import show as rio_show
from matplotlib.pyplot import subplots, savefig
from numpy import zeros
from matplotlib.colors import LinearSegmentedColormap
from skimage.draw import line
from numpy import column_stack

# open the elevation data file
with rio_open("../data/helvellyn/Helvellyn-50.asc") as d:

    # read the data out of band 1 in the dataset
    dem_data = d.read(1)
    # convert coordinates to image space    
    row, col = d.index(334241, 515107)
    print(line(row, col, row, col+10))
    print(column_stack(line(row, col, row, col+10)))

  
    # plot the dataset
    fig, my_ax = subplots(1, 1, figsize=(16, 10))
    
    # add the DEM
    rio_show(
      dem_data,
      ax=my_ax,
      transform = d.transform,
    )
    
    # create a new 'band' of raster data the same size
    output = zeros((d.height, d.width))
    output[row][col] = 1
    
    for row, col in column_stack(line(row, col, row, col+10)):
        # colour in red
        output[row][col] = 1
    
    # add the drawing layer
    rio_show(
        output,
        cmap = LinearSegmentedColormap.from_list('binary', [(0, 0, 0, 0), (1, 0, 0, 0.5)], N=2),
        ax=my_ax,
        transform=d.transform   
        )
    
    
    # save the resulting map
    savefig('./out/bresenham.png', bbox_inches='tight')