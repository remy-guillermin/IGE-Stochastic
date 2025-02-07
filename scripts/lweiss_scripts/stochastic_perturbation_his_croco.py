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
from scipy.stats import gaussian_kde

#data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
data = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

member = '001'
simu = 'run_swio_stoens_2017_02'

### grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk == 0, msk, np.nan)
area = 1 / (g['pm'] * g['pn'])
area = area.where(msk == 1, np.nan)
#sys.exit()
g.close()

# Définir la période pour calculer la moyenne temporelle
start_time = '2017-01-31T12:00:00' # '2019-01-01T12:00:00'
end_time = '2017-08-22T12:00:00'

# Remplacer les valeurs hors domaine par NaN
ds1 = xr.open_dataset(data + simu + f'/{member}swiose_his.nc')  # , engine='h5netcdf')
p_sto1 = ds1['tmpout'][:, -1, :, :]
p_sto1 = p_sto1.where((msk == 1) & (p_sto1 != 9.96921e+36), np.nan)
area2 = xr.DataArray(np.tile(area.expand_dims('time').values, (ds1.dims['time'], 1, 1)),
                    dims=['time', 'eta_rho', 'xi_rho'],
                    coords={'time': ds1.time, 'eta_rho': area.eta_rho, 'xi_rho': area.xi_rho})
p_sto_t_mean1 = np.nansum(p_sto1, axis=(1, 2)) / np.nansum(area2[:,:,:].values, axis=(1, 2))
p_sto_t_mean1 = p_sto1[:,171,83]
# Sélectionner la période
p_sto_period1 = p_sto1.sel(time=slice(start_time, end_time))
# Calculer la moyenne temporelle sur la période sélectionnée
p_sto_mean1 = p_sto_period1.mean(dim='time', skipna=True)

# ds2 = xr.open_dataset(data + simu + '_0002/002swiose_his.nc')  # , engine='h5netcdf')
# p_sto2 = ds2['tmpout'][:, -1, :, :]
# p_sto2 = p_sto2.where(msk == 1, np.nan)
# Sélectionner la période
# p_sto_period2 = p_sto2.sel(time=slice(start_time, end_time))
# Calculer la moyenne temporelle sur la période sélectionnée
# p_sto_mean2 = p_sto_period2.mean(dim='time', skipna=True)
# p_sto_bias = p_sto_mean2 - p_sto_mean1

################################################################################################################
fig = plt.figure(figsize=(10, 5)) #, constrained_layout=True)

ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
cmap = cmcrameri.cm.vik
a = 0.8 # -0.1
b = 1.2 # 0.1
c = 9 # 10 # int((b - a) + 1)
levels = np.linspace(a, b, c)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)
# Map
pcm = ax1.pcolormesh(lon, lat, p_sto_mean1, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax1.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax1.contourf(lon, lat, msk_inv, colors='lightgray')
space_mean=p_sto_mean1.mean(skipna=True)
ax1.text(0.02, 0.98, f'Diffusive Operator\nwith 100 passes\nmean={space_mean:.2f}', transform=ax1.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3",edgecolor='black',facecolor='white',alpha=0.5))
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}
ax1.set_title(f"Spatial Correlation - member{member} - {start_time[:13]}", size=9)
### colorbar
cb = fig.colorbar(pcm, ax=ax1, orientation='vertical', fraction=0.046, pad=0.04, label='')
colorbaryticks = np.linspace(a, b, int(c/2)+1) # 5  # levels[::2]
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),2), fontsize=8)

# Subplot 2: Moyenne temporelle des valeurs spatiales (col2, ligne1)
ax2 = fig.add_subplot(2, 2, 2)
# Tracer l'évolution temporelle de la moyenne spatiale
ax2.plot(np.arange(0, len(ds1.time.data)), p_sto_t_mean1, color='k',linewidth=0.8,linestyle='-')
ax2.set_title("Time Correlation", fontsize=10)
# ax2.set_xlabel("Time", fontsize=9)
ax2.set_xticks(np.arange(0, len(ds1.time.data), 24*30))
ax2.set_xticklabels(np.arange(1, int(len(ds1.time.data)/24) + 2, 30), rotation=0)
ax2.grid(True, linestyle='--', alpha=0.3)
ax2.set_xlim(0, len(ds1.time.data))
ax2.set_ylim(0, 2)
# ax2.set_ylabel("Mean value", fontsize=9)
time_mean=p_sto_t_mean1[1:].mean(skipna=True)
ax2.text(0.02, 0.96, f'Autoregressive Process order 2\n30min time interp.\ntime correlation = 10 days\nmean={time_mean:.2f}', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3",edgecolor='black',facecolor='white',alpha=0.5))

# Subplot 3: Histogramme avec la PDF des valeurs spatiales (col2, ligne2)
ax3 = fig.add_subplot(2, 2, 4)
# Aplatir les données spatiales sur tous les pas de temps pour construire un histogramme
flattened_data = p_sto1[1:,:,:].values.flatten()
# Tracer l'histogramme et la PDF
ax3.hist(flattened_data, bins=1000, density=True, color='b', alpha=0.6)
ax3.set_title("Marginal distribution", fontsize=10)
# ax3.set_xlabel("stochastic parameter value", fontsize=9)
ax3.set_ylabel("PDF", fontsize=9)
ax3.grid(True, linestyle='--', alpha=0.3)
# kde_nan = flattened_data[np.isfinite(flattened_data)]
# kde = gaussian_kde(kde_nan)
# x_vals = np.linspace(np.min(kde_nan), np.max(kde_nan), 1000)  # Valeurs x pour la PDF
# pdf_vals = kde(x_vals)
# ax3.plot(x_vals, pdf_vals, color='b', lw=2) #, label='PDF')
dist_mean=np.nanmean(p_sto1[1:,:,:].values)
dist_min=np.nanmin(p_sto1[1:,:,:].values)
dist_max=np.nanmax(p_sto1[1:,:,:].values)
dist_std=np.nanstd(p_sto1[1:-1,:,:].values)
dist_pts=len(flattened_data)
ax3.text(0.8, 0.96, f'Lognormal\nmean={dist_mean:.2f}\nstd={dist_std:.2f}\nmin={dist_min:.2f}\nmax={dist_max:.2f}', transform=ax3.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3",edgecolor='black',facecolor='white',alpha=0.5))
ax3.set_xlim(0, dist_max)

### SAVE
plt.tight_layout()
plt.savefig(path_fig + simu + f'/StoDist_{simu[:-8]}{member}_{end_time[:13]}.png', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()
ds1.close()
# ds2.close()