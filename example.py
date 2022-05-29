# %%
from connectivity import Sector, Connectivity
import netCDF4 as nc

# %%
# get trajectories data
ncfiles = r"/mnt/c/users/bassi/krill_transport_wap_palmer_canyon_shelf_F2_20200101_30_leeway_50.nc"
ds = nc.Dataset(ncfiles)
lons = ds['lon'][:]
lats = ds['lat'][:]

# %%
## Example for creating from shapefile ###
# create sectors from shapefile
shapefile = r"/mnt/c/users/bassi/Downloads/wap_chlora_202001.shp"
shapefile_sectors = Sector.sectors_from_shapefile(shapefile)

# %%
shapefile_model = Connectivity(lons, lats, shapefile_sectors)

# %%
shapefile_model.plot(shapefile_model.n_settle_values, title="Nsettle")

# %%
### Example for creating from reactangular domain ###
# create sectors from domain
domain = [-60, -67, -62, -6]
domain_sectors = Sector.create_sectors(domain, 0.1)

# %%
domain_model = Connectivity(lons, lats, domain_sectors)
# %%
domain_model.plot(domain_model.n_settle_values, squares=True, delta=0.1, corners=domain, title='Nsettle')
# %%
domain_model.plot(domain_model.n_supply_values, squares=True, delta=0.1, corners=domain)

