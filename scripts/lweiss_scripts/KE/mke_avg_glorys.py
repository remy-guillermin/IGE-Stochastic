## LIB
import numpy as np
import pandas as pd
import netCDF4 
import xarray as xr
import matplotlib as mpl 
import matplotlib.pyplot as plt
import cmocean 
import cmcrameri
import sys
import cartopy.crs as ccrs
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.patches as patches
##

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/DATA/GLORYS/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'cmems_glorys'

### open netcdf file an variables
g = xr.open_dataset('/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIO/run_swio_ref_2017_2023_restart1/swiose_grid.nc')
lon_g = g['lon_rho'][:, :]
lat_g = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk==0, msk, np.nan)
g.close()

ds = xr.open_dataset(data + '/raw_cli_mercator_Y2019.nc')
lon = ds['longitude'][:]
lat = ds['latitude'][:]
u = ds['uo'][:, 0, :, :]  # Sélectionne le premier niveau s_rho
v = ds['vo'][:, 0, :, :]

# Remplacer les valeurs hors domaine par NaN
#fill_value = 9.96921e+36
#u = u.where((u != fill_value), np.nan)
#v = v.where((v != fill_value), np.nan)

# Définir la période pour calculer la moyenne temporelle
start_time = '2019-01-01' # '2019-01-01T12:00:00'
end_time = '2019-12-31'

# Sélectionner la période
# u_period = u.sel(time=slice(start_time, end_time))
# v_period = v.sel(time=slice(start_time, end_time))

u_mean = u.mean(dim='time')
v_mean = v.mean(dim='time')

# ds.close()

uu = units.Quantity(u_mean[:,:].data, "m/s")
vv = units.Quantity(v_mean[:,:].data, "m/s")

mke_mean = 0.5 * (uu ** 2 + vv ** 2)
# mke_mean = np.nanmean(mke, axis=0)

print('U, V, mke, lon, lat, msk \n', mke_mean.shape, lon.shape, lat.shape, msk.shape)

### plot
boxes = [(52,60,-24,-16), (41,47,-15,-8), (48,60,-4,3), (36.5,42.5,-28,-19)]
names = ['Mascarene','Mayotte-Comores','Equator','South Moz']

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
cmap = cmcrameri.cm.roma_r
# a = 0
# b = 0.5
# c = 11
# levels = np.linspace(a, b, c*2-1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.LogNorm(vmin=1e-2, vmax=1e0)
# norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, mke_mean, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon_g, lat_g, msk, colors='k', linewidths=0.1)
ax.contourf(lon_g, lat_g, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}
#ax.coastlines(resolution='10m')
ax.set_title(f"MKE GLORYS {start_time} to {end_time}", size=10)

# Ajouter les boîtes en contours fins noirs
# for i, box in enumerate(boxes):
#     x_min, x_max, y_min, y_max = box
#     width = x_max - x_min
#     height = y_max - y_min
#     rect = patches.Rectangle((x_min, y_min), width, height, linewidth=0.7, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree())
#     ax.add_patch(rect)
# 
#     # Ajouter le nom au centre de chaque boîte
#     ax.text(x_min + 0.2, y_max - 1.2, names[i],
#             # horizontalalignment='center', verticalalignment='center',
#            transform=ccrs.PlateCarree(), fontsize=6, color='black') #, weight='bold')


### colorbar
cb = fig.colorbar(pcm, ax=ax, label='Mean Kinetic Energy $[m^{2} s^{-2}]$')
# colorbaryticks = np.linspace(a, b, c)  # levels[::2]
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
# cb.set_ticks(colorbaryticks)
# cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),1), fontsize=8)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=8)
text.set_font_properties(font)

### SAVE
plt.savefig(path_fig + simu + f'/mke_mean_{simu[:13]}_{start_time}_{end_time}.png', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()
