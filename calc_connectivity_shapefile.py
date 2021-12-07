import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift, timedelta
from calc_connectivity import *

# %%
# load data
ncfiles = r"/mnt/c/users/bassi/data 2 - Copy/*.nc"
shapefile = "/mnt/c/users/bassi/research/gis/איזור_חיפושי_שמורות_ימיות/איזורי_חיפוש_שמורות_ימיות.shp"

# %%
# create reader and initialize objects
reader = reader_netCDF_CF_generic.Reader(ncfiles)
o = OceanDrift(loglevel=30)
o.set_config('seed:ocean_only', False)
o.add_reader([reader])
domain = [33.6, 34.5, 31.5, 34]
sectors = sectors_from_shapefile(shapefile)

# %%
# prepare and run simulation
first_lonlats, last_lonlats, supply, sectors = run_for_connectivity_polygons(o,
                                                                             duration=timedelta(days=5),
                                                                             time_step=timedelta(hours=12),
                                                                             d_supply=500,
                                                                             sectors=sectors,
                                                                             time=[reader.start_time + timedelta(0),
                                                                                   reader.start_time + timedelta(
                                                                                       4) + timedelta(0)], )

# %%
# get connectivity matrix
C = calculate_connectivity_matrix(first_lonlats, last_lonlats, supply_org=supply,
                                  sectors=sectors, polygons=True)

# %%
# get connectivity values
Nsettle_values, Nsupply_values, R_values, Srec_values, Ssup_values, Sector_categories = categorize_sectors(
    sectors=sectors, C=C, cat_threshold=4)

# %%
# plot connectivity matrix
import matplotlib.pyplot as plt
import numpy as np

cmap = plt.get_cmap('viridis')
cmap.set_under(color='white')
plt.imshow(C, cmap=cmap, interpolation='nearest', vmin=0.001)
plt.show()

# plot connectivity categories
plot_connectivity_polygons(sectors=sectors, sector_categories=Sector_categories, title="Categories", vmax=8, lut=8)
