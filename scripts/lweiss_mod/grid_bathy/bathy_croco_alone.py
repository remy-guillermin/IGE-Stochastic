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
from scipy.interpolate import interp1d
##

work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/Grid/'
grid = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'

#######################################################################################
### open grid
ds = xr.open_dataset(grid)
h = ds['h'][:, :]
lon = ds['lon_rho'][:, :]
lat = ds['lat_rho'][:, :]
msk = ds['mask_rho'][:, :]
#sys.exit()
ds.close()

msk_inv = np.where(msk==0, msk, np.nan)

#######################################################################################
### plot
#######################################################################################
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
cmap = cmocean.cm.deep_r
a = 0
b = 5000
c = int((b-a)*2/1000+1)
levels = np.linspace(a, b, c*2-1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, h, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cs = ax.contour(lon, lat, h, levels=np.arange(0, 6000, 1000), colors='k',
                linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs, fmt='%d', inline=True, fontsize=5)
cs2 = ax.contour(lon, lat, h, levels=np.arange(-300, 600, 400), colors='red',
                linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs2, fmt='%d', inline=True, fontsize=5)
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

#### zoom ####
# pcm = ax.pcolormesh(lon[zoom], lat[zoom], h[zoom[0], zoom[1]], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
# cs = ax.contour(lon[zoom], lat[zoom], h[zoom[0], zoom[1]], levels=np.arange(0, 5000, 500), colors='k',
#                linestyles='dashed', linewidths=0.2, transform=ccrs.PlateCarree())
# plt.clabel(cs, fmt='%d', inline=True, fontsize=4)
# cs2 = ax.contour(lon[zoom], lat[zoom], h[zoom[0], zoom[1]], levels=np.arange(0, 500, 100), colors='red',
#                linestyles='dashed', linewidths=0.2, transform=ccrs.PlateCarree())
# plt.clabel(cs2, fmt='%d', inline=True, fontsize=4)
# ax.contour(lon[zoom], lat[zoom], msk[zoom], colors='k', linewidths=0.1)
# ax.contourf(lon[zoom], lat[zoom], msk_inv[zoom], colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 6, 'color': 'k'}
gl.ylabel_style = {'size': 6, 'color': 'k'}
#ax.coastlines(resolution='10m')
#ax.set_title('Sea Surface Temperature')

# Add vertical and horizontal grid lines (based on lon, lat)
step = 10
for xi in range(0, lon.shape[1], step):  # Loop through the xi (longitude) points
    ax.plot(lon[:, xi], lat[:, xi], color='gray', linewidth=0.3, alpha=0.6, transform=ccrs.PlateCarree())  # Vertical lines
for eta in range(0, lat.shape[0], step):  # Loop through the eta (latitude) points
    ax.plot(lon[eta, :], lat[eta, :], color='gray', linewidth=0.3, alpha=0.6, transform=ccrs.PlateCarree())  # Horizontal lines

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='Bathymetry (m)')
colorbaryticks = np.linspace(a, b, c)  # levels[::2]
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=6)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=6)
text.set_font_properties(font)

### SAVE
plt.show()

