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

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

g = xr.open_dataset(grid)
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

d = xr.open_dataset(data + simu + 'swio_avg_2017.nc')
zeta = d['zeta'][:, :, :]
d.close()

start_time = '2017-01-31'
end_time = '2017-01-31'

zeta_period = zeta.sel(time=slice(start_time, end_time))
print('Height: ', zeta_period.shape)

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
zeta = zeta.where((zeta != fill_value), np.nan)

# Moyenne sur start - end
zeta_mean = zeta_period.mean(dim='time')

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_title(f"Mean SSH SWIO from {start_time} to {end_time}", size=9)

cmap = cmcrameri.cm.roma_r
a = 0
b = 1
c = 11 # int((b - a) + 1)
levels = np.linspace(a, b, c * 2 -1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, zeta_mean, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='Sea Surface Height [m]')
colorbaryticks = np.linspace(a, b, c) 
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height]) 
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),1), fontsize=8)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=8)
text.set_font_properties(font)

plt.show()