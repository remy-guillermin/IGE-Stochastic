# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import glob
import os
from matplotlib.gridspec import GridSpec

work = '/lus/work/CT1/c1601279/rguillermin'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'

stoens = 'run_swio2_stoens30_gls'

grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

lat_index = [430, 280]
SWIO = (25, 69, -36, 7)

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
mask_rho = g.mask_rho
g.close()
print('Grid loaded')

ds = xr.open_dataset(os.path.join(work, 'MLD', stoens, '001swiose_mld.nc'))

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

date = ds.isel(time=-1).time.data.astype('datetime64[D]')

fig.suptitle(f"Mixing density layer - {date}")

ax.set_extent(SWIO)   
ax.coastlines(resolution='50m')
ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'k'}
gl.ylabel_style = {'size': 10, 'color': 'k'}

data = ds.isel(time=-1).mld

pcm = ax.pcolormesh(lon, lat, data, cmap=cmocean.cm.deep, transform=ccrs.PlateCarree(), vmin=-110, vmax=0)
cb = fig.colorbar(pcm, ax=ax, label='Mixing layer depth [m]')

fig.tight_layout()
fig.savefig(f'{figures}/mld_interp_{date}.png')
plt.close()
