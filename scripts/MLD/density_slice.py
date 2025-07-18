# LIBRAIRIES
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cmocean
import cartopy.crs as ccrs
import glob
import os
from matplotlib.gridspec import GridSpec
%matplotlib inline



figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
data = '/lus/store/CT1/c1601279/lweiss/run_croco/SWIO/run_swio2_deter2_2017_2023'
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

CRIT = 0.03
lat_index = 280 #430



g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
mask_rho = g.mask_rho
g.close()

h = h.where(mask_rho, 0)

depth_level = xr.open_dataset(depth_lvl)



ds = xr.open_dataset(os.path.join(data, 'swio_avg_2019.nc'))[['rho']].isel(time=100)
ds = ds.where(mask_rho, np.nan)
ds



np.tile(lon.isel(eta_rho=lat_index).values, (50, 1)).shape, depth_level.isel(eta_rho=lat_index).depth_level.data.shape



fig, ax = plt.subplots(1, 1, figsize=(10,4))

# fig.suptitle('Zonal slice of the density field')

land_color = ccrs.cartopy.feature.COLORS['land']

ax.fill_between(lon[lat_index, :], -h[lat_index, :].values, y2=min(-h[lat_index, :]), color=land_color, label='Land')
ax.plot(lon[lat_index, :], -h[lat_index, :].values, color='black', linewidth=0.5)

rho = ds.rho.isel(eta_rho=lat_index) + 1025

pcm = ax.pcolormesh(np.tile(lon.isel(eta_rho=lat_index).values, (50, 1)), depth_level.isel(eta_rho=lat_index).depth_level.data, rho.data, cmap = cmocean.cm.dense, shading='auto')

cb = fig.colorbar(pcm, ax=ax, orientation='vertical', pad=0.02)
cb.set_label('Density [kg/m³]')

ax.text(49.3, 10, 'Madagascar', ha='center', va='bottom')
ax.text(43.7, 10, 'Comoros', ha='center', va='bottom')
    
ax.set_xlim(40, 65)
ax.set_ylim(bottom=-500, top=0)

ax.set_xlabel('Longitude')
ticks = ax.get_xticks()
ax.set_xticks(ticks)
ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])

ax.set_ylabel('Depth')
ticks = ax.get_yticks()
ax.set_yticks(ticks)
ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])

ax.legend(loc='upper right', framealpha=1)

fig.tight_layout()
fig.savefig(os.path.join(figures, 'density_slice.png'), dpi=300, transparent=True)
plt.show()