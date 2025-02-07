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
from scipy.stats import pearsonr

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_stogen_1mth'
simu = 'run_swio_stogen_1yr_diff_100p_10jo1ud12'

lon_index_a, lat_index_a = 83, 171
lon_index_b, lat_index_b = 290, 151
lon_index_c, lat_index_c = 231, 328

### open netcdf file and variables
### open grid
g = xr.open_dataset(data + simu + '_0001/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk == 0, msk, np.nan)
#sys.exit()
lon_a, lat_a = g['lon_rho'][lat_index_a, lon_index_a], g['lat_rho'][lat_index_a, lon_index_a]
lon_b, lat_b = g['lon_rho'][lat_index_b, lon_index_b], g['lat_rho'][lat_index_b, lon_index_b]
lon_c, lat_c = g['lon_rho'][lat_index_c, lon_index_c], g['lat_rho'][lat_index_c, lon_index_c]
g.close()

lon_index_a = 83
lat_index_a = 171

lon_index_b = 290
lat_index_b = 151

lon_index_c = 231
lat_index_c = 328

# Définir la période pour calculer la moyenne temporelle
start_time = '2017-09-01T12:00:00' # '2019-01-01T12:00:00'
end_time = '2017-09-01T12:00:00'

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36

ds1 = xr.open_dataset(data + simu + '_0001/001swiose_his.nc')  # , engine='h5netcdf')
zeta1 = ds1['zeta'][:, :, :]
zeta1 = zeta1.where(zeta1 != fill_value, np.nan)
# Sélectionner la période
zeta_period1 = zeta1.sel(time=slice(start_time, end_time))
# Calculer la moyenne temporelle sur la période sélectionnée
zeta_mean1 = zeta_period1.mean(dim='time', skipna=True)

ds2 = xr.open_dataset(data + simu + '_0002/002swiose_his.nc')  # , engine='h5netcdf')
zeta2 = ds2['zeta'][:, :, :]
zeta2 = zeta2.where(zeta2 != fill_value, np.nan)
# Sélectionner la période
zeta_period2 = zeta2.sel(time=slice(start_time, end_time))
# Calculer la moyenne temporelle sur la période sélectionnée
zeta_mean2 = zeta_period2.mean(dim='time', skipna=True)

# Tracer la figure représentant la moyenne temporelle
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

cmap = cmcrameri.cm.roma_r #roma_r #vik
a = 0 #-0.2 # -0.1
b = 1.2 #1.2 # 0.1
c = 13 #15 # 10 # int((b - a) + 1)
levels = np.linspace(a, b, c)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, zeta_mean2, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

# Dessiner les points P1 et P2
ax.plot(lon_a, lat_a, marker='o', color='k', markersize=1, transform=ccrs.PlateCarree())
ax.plot(lon_b, lat_b, marker='o', color='k', markersize=1, transform=ccrs.PlateCarree())
ax.plot(lon_c, lat_c, marker='o', color='k', markersize=1, transform=ccrs.PlateCarree())

# Ajouter les étiquettes P1 et P2
# ax.text(lon_a + 0.5, lat_a, 'P1', color='k', transform=ccrs.PlateCarree(), fontsize=9) #, weight='bold')
# ax.text(lon_b + 0.5, lat_b, 'P2', color='k', transform=ccrs.PlateCarree(), fontsize=9) #, weight='bold')
# ax.text(lon_c + 0.5, lat_c, 'P3', color='k', transform=ccrs.PlateCarree(), fontsize=9) #, weight='bold')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}
# ax.coastlines(resolution='10m')
ax.set_title(f"SSH - member 2 - {start_time[:13]}", size=9)

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='Sea Surface Height [m]')
colorbaryticks = np.linspace(a, b, 7) #int(c/2)+1) # 5  # levels[::2]
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),2), fontsize=8)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=8)
text.set_font_properties(font)

### SAVE
plt.savefig(path_fig + simu + f'/ssh_{simu}_002_{start_time[:13]}.png', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()

ds1.close()
ds2.close()
