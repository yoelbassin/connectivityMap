import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
from opendrift.models.oceandrift import OceanDrift, timedelta
from opendrift.models.plastdrift import PlastDrift
from os import listdir
from calc_connectivity import *

# %% load data
ncfiles = r"/mnt/c/users/bassi/data 2 - Copy/*.nc"
# for lazy:
# root = r"/mnt/c/users/bassi/data 2 - Copy/"
# ncfiles_separate=[root + f for f in listdir(root) if f.split('.')[-1]=='nc']


# %%
# get the ports coordinates for the simulation
# source:
# https://geonode.wfp.org/layers/esri_gn:geonode:wld_trs_ports_wfp
# download:
# https://geonode.wfp.org/geoserver/wfs?typename=geonode%3Awld_trs_ports_wfp&outputFormat=csv&version=1.0.0&request=GetFeature&service=WFS
ports_coords = []
# get ports locations
from csv import reader
# skip first line i.e. read header first and then iterate over each row od csv as a list
with open(r"./data/wld_trs_ports_wfp.csv", 'r') as read_obj:
    csv_reader = reader(read_obj)
    header = next(csv_reader)
    # Check file as empty
    if header is not None:
        # Iterate over each row after the header in the csv
        for row in csv_reader:
            # row variable is a list that represents a row in csv
            ports_coords.append((float(row[12]), float(row[11])))

# %%
# create reader and initialize objects
reader_landmask = reader_global_landmask.Reader()
ncreader = reader_netCDF_CF_generic.Reader(ncfiles)
o = PlastDrift(loglevel=30)
o.set_config('seed:ocean_only', False)
# for lazy:
# o.add_reader_from_list(ncfiles_separate, lazy=True)
o.add_reader([ncreader])
o.add_reader([reader_landmask])
domain = [20, 37, 30, 41]
delta = 0.2

# %%
# prepare and run simulation
first_lonlats, last_lonlats, supply, sectors = run_for_connectivity_domain(o,
                                                                           duration=timedelta(days=35),
                                                                           time_step=timedelta(hours=24),
                                                                           domain=domain,
                                                                           supply_coords=ports_coords,
                                                                           delta=delta,
                                                                           d_supply=3500,
                                                                           time=[ncreader.start_time + timedelta(0),
                                                                                 ncreader.start_time + timedelta(
                                                                                     34) + timedelta(0)],
                                                                           _radius=500,
                                                                           time_step_output=timedelta(days=2))

# %%
# get connectivity matrix
C = calculate_connectivity_matrix(first_lonlats, last_lonlats, supply_org=supply, sectors=sectors)

# %%
# get connectivity values
Nsettle_values, Nsupply_values, R_values, Srec_values, Ssup_values, Sector_categories = categorize_sectors(
    sectors=sectors, C=C)

# %%
# plot Nsettle
plot_connectivity_map(sectors=sectors, sector_categories=Nsettle_values, title="Nsettle")

# %%
# plot Nsupply
plot_connectivity_map(sectors=sectors, sector_categories=Nsupply_values, title="Nsupply")

# %%
# plot R values
plot_connectivity_map(sectors=sectors,sector_categories=R_values, title="Settling/Supply factor")

# %%
# plot connectivity categories
plot_connectivity_map(sectors=sectors, sector_categories=Sector_categories, title="Categories", vmax=8, lut=8)
