# %%
from connectivity import Sector, Connectivity
import netCDF4 as nc

domain = [-60, -67, -62, -6]

sectors = Sector.create_sectors(domain, 0.1)

ncfiles = r"/mnt/c/users/bassi/krill_transport_wap_palmer_canyon_shelf_F2_20200101_30_leeway_50.nc"

ds = nc.Dataset(ncfiles)

lons = ds['lon'][:]
lats = ds['lat'][:]
# %%
c = Connectivity(lons, lats, sectors)

# %%
c.plot(c.n_settle_values, squares=True, delta=0.1, corners=domain, title='Nsettle')
# %%
c.plot(c.n_supply_values, squares=True, delta=0.1, corners=domain)