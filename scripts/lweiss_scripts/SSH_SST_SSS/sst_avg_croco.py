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

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_0003_2017'

### open netcdf file and variables
### open grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
#sys.exit()
g.close()

ds = xr.open_dataset(data + simu + '/swiose_avg.nc')  # , engine='h5netcdf')
temp = ds['temp'][:, -1, :, :]  # Sélectionne le premier niveau s_rho
# lon = ds['lon_rho'][:, :]
# lat = ds['lat_rho'][:, :]
# msk = ds['mask_rho'][:, :]

msk_inv = np.where(msk == 0, msk, np.nan)

# Print the shape for debugging purposes
print(temp.shape, lon.shape, lat.shape, msk.shape)

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
temp = temp.where((temp != fill_value) & (ds.mask_rho == 1), np.nan)

# Définir la période pour calculer la moyenne temporelle
start_time = '2017-10-01' # '2019-01-01T12:00:00'
end_time = '2017-12-31'

# Sélectionner la période
temp_period = temp.sel(time=slice(start_time, end_time))

# Calculer la moyenne temporelle sur la période sélectionnée
temp_mean = temp_period.mean(dim='time', skipna=True)

# Sauvegarder temp_mean comme NetCDF
temp_mean.to_netcdf(path_fig + simu + f'/Mean_SST/sst_mean_{simu}_{start_time}_{end_time}.nc')

sys.exit()

# Tracer la figure représentant la moyenne temporelle
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

cmap = cmocean.cm.thermal
a = 20
b = 30
c = int((b - a) + 1)
levels = np.linspace(a, b, c * 2 - 1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, temp_mean, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}
# ax.coastlines(resolution='10m')
ax.set_title(f"Mean SST SWIO from {start_time} to {end_time}", size=9)

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='Sea Surface Temperature °C')
colorbaryticks = np.linspace(a, b, c)  # levels[::2]
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=8)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=8)
text.set_font_properties(font)

### SAVE
plt.savefig(path_fig + simu + f'/sst_mean_{simu}_{start_time}_{end_time}.png', dpi=300, bbox_inches='tight')
plt.close()

ds.close()

