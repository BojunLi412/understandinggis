# -*- coding: utf-8 -*-
"""
Created on Wed May  5 10:04:54 2021

@author: 81465
"""

from geopandas import read_file
from sys import exit
from fiona.errors import DriverError

with rio_open('../data/natural-earth/ne_10m_admin_0_countries.shp') as dem: