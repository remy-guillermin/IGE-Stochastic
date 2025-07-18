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

h = h.where(h != 50, 0)

depth_level = xr.open_dataset(depth_lvl)

fig, axs = plt.subplots(1, 2, figsize=(12,4), gridspec_kw={'width_ratios': [2, 1]})

levels = depth_level.isel(eta_rho=lat_index)

for ax in axs:

    land_color = ccrs.cartopy.feature.COLORS['land']
    
    for level in levels.depth_level.values:
        ax.plot(lon[lat_index, :], level, color='gray', linewidth=0.5)

    ax.fill_between(lon[lat_index, :], -h[lat_index, :].values, y2=min(min(-h[lat_index, :]), -5000), color=land_color, label='Land', zorder=2)
    ax.plot(lon[lat_index, :], -h[lat_index, :].values, color='black', linewidth=0.5)
    
ax = axs[0]
ax.text(49.3, 10, 'Madagascar', ha='center', va='bottom')
ax.text(43.7, 10, 'Comoros', ha='center', va='bottom')

ax.set_xlim(40, 65)
ax.set_ylim(-5000, 0)

ax.legend(loc='lower left', framealpha=1)

ax = axs[1]
ax.set_xlim(47, 52)
ax.set_ylim(-1000, 0)

for ax in axs:
    ax.set_xlabel('Longitude')
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])

    ax.set_ylabel('Depth')
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])

fig.tight_layout()
fig.savefig(os.path.join(figures, 'levels_slice.png'), dpi=300)
plt.show()
