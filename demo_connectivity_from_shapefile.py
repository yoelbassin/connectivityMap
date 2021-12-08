# This file is part of ConnectivityMap module
#
# ConnectivityMap is a free module for creating a connectivity map 
# and matrix from a lagrangian trajectory model
#
# Copyright 2021, Erick Fredj and Yoel Bassin, JCT, Israel


import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift, timedelta
from calculate_connectivity import *
"""
This demo is an example for the use of ConnectivityMap using a shape file.
The demo receives a shapefile with polygons that describe sea sectors, and using
a lagrangian trajectory model it generates a connectivity matrix and map.
"""

# %%
# load data
ncfiles = r"/mnt/c/users/bassi/data 2 - Copy/*.nc"
shapefile = "/mnt/c/users/bassi/research/gis/איזור_חיפושי_שמורות_ימיות/איזורי_חיפוש_שמורות_ימיות.shp"

# %%
# create reader and initialize objects
reader = reader_netCDF_CF_generic.Reader(ncfiles)
o = OceanDrift(loglevel=30) # Passive opendrift elements
o.set_config('seed:ocean_only', False)
o.add_reader([reader])
domain = [33.6, 34.5, 31.5, 34]
sectors = sectors_from_shapefile(shapefile)
supply_amount = 3000

# %%
# prepare and run simulation
first_lonlats, last_lonlats, supply, sectors = run_for_connectivity_polygons(o,
                                                                             duration=timedelta(days=30),
                                                                             time_step=timedelta(hours=12),
                                                                             supply_ammount=supply_amount,
                                                                             sectors=sectors,
                                                                             time=[reader.start_time + timedelta(0),
                                                                                   reader.start_time + timedelta(
                                                                                       29) + timedelta(0)], )

# %%
# get connectivity matrix
C = calculate_connectivity_matrix(first_lonlats, last_lonlats, supply=supply,
                                  sectors=sectors, polygons=True)

# %%
# get connectivity values
Nsettle_values, Nsupply_values, R_values, Srec_values, Ssup_values, Sector_categories = categorize_sectors(
    sectors=sectors, C=C, cat_threshold=4)

# %%
# plot connectivity matrix
import matplotlib.pyplot as plt
import numpy as np

plot_connectivity_network(sectors=sectors, connectivity_matrix=C/supply_amount,
 sector_categories=Sector_categories)

cmap = plt.get_cmap('viridis')
cmap.set_under(color='white')
fig, ax = plt.subplots(1,1)

img = ax.imshow(C/supply_amount, cmap=cmap, interpolation='nearest', vmin=0.001)

fig.colorbar(img)
plt.show()

# plot connectivity categories
plot_connectivity_polygons(sectors=sectors, sector_categories=Sector_categories, title="Categories", vmax=8, lut=8, numbering=True)
