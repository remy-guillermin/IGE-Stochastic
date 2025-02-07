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

scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'
grid = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swiose.nc'
simu = 'run_swiose_0004'

boxes = [(35,40,-23,-19)] # MoC

#######################################################################################
### open grid
ds = xr.open_dataset(grid)
h = ds['h'][:, :]
lon = ds['lon_rho'][:, :]
lat = ds['lat_rho'][:, :]
msk = ds['mask_rho'][:, :]
#sys.exit()
ds.close()

### open netcdf file an variables
ds = xr.open_dataset(scratch + simu + '/swiose_his.nc') # , engine='h5netcdf')
temp = ds['temp'][:, :, :, :]
#h = ds['h'][:, :]  # bathymetry at RHO-points
#lon = ds['lon_rho'][:, :]
#lat = ds['lat_rho'][:, :]
#msk = ds['mask_rho'][:, :]
s_rho = ds['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = ds['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = ds['hc'].values # S-coordinate parameter, critical depth

# sys.exit()
ds.close()

msk_inv = np.where(msk==0, msk, np.nan)

#######################################################################################
# Calcul z0 a nonlinear vertical transformation 
def calc_depth(s, Cs, hc, h):
    N = len(s_rho)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h ## (hc * s[k] + h * Cs[k])
    return depth

# Calcul de la profondeur
depth_sigma = calc_depth(s_rho, Cs_rho, hc, h)

#######################################################################################
#### h transect
#######################################################################################
lat_index = 243 # 249
lon_index = 162 # 195

# along longitude axis
fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})

ax[0].plot(lat[:, lon_index].values, -h[:, lon_index].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[0].fill_between(lat[:, lon_index], -h[:, lon_index].values, y2=min(-h[:, lon_index]), color='lightgrey')

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[0].plot(lat[:, lon_index], depth_sigma[k, :, lon_index], color='grey', linestyle='-', linewidth=0.5)
ax[0].set_xlim(-16.5, 3)
ax[0].set_ylim(np.min(-h[:, lon_index].values), 5)
ax[0].set_xlabel('Latitudes along ' + str(np.round(lon[lat_index, lon_index].values,2)) + '°E Longitude')
ax[0].set_ylabel('Depth h (m)')
ax[0].grid(linestyle='--',linewidth=0.3)

### subplot
ax[1].plot(lat[:, lon_index].values, -h[:, lon_index].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[1].fill_between(lat[:, lon_index], -h[:, lon_index].values, y2=min(-h[:, lon_index]), color='lightgrey')

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[1].plot(lat[:, lon_index], depth_sigma[k, :, lon_index], color='grey', linestyle='-', linewidth=0.5)
ax[1].set_xlim(-13, -12.7)
ax[1].set_ylim(-100, 1)
ax[1].grid(linestyle='--',linewidth=0.3)

plt.savefig(path_fig + 'transect_zoom_h_lon_' + str(np.round(lon[lat_index, lon_index].values,2)) + '_' + simu + '.png', dpi=300, bbox_inches='tight')
# plt.show()

#######################################################################################
# along latitude axis
fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})
ax[0].plot(lon[lat_index,:].values, -h[lat_index,:].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[0].fill_between(lon[lat_index,:], -h[lat_index,:].values, y2=min(-h[lat_index, :]), color='lightgrey')
# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[0].plot(lon[lat_index, :], depth_sigma[k, lat_index, :], color='grey', linestyle='-', linewidth=0.5)
ax[0].set_xlim(40, 49)
ax[0].set_ylim(-4000, 0)
ax[0].set_xlabel('Longitudes along ' + str(np.round(lat[lat_index, lon_index].values,2)) + '°S Latitude')
ax[0].set_ylabel('Depth (m)')
ax[0].grid(linestyle='--',linewidth=0.3)

### subplot
ax[1].plot(lon[lat_index,:].values, -h[lat_index,:].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[1].fill_between(lon[lat_index,:], -h[lat_index,:].values, y2=min(-h[lat_index, :]), color='lightgrey')
# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[1].plot(lon[lat_index, :], depth_sigma[k, lat_index, :], color='grey', linestyle='-', linewidth=0.5)
ax[1].set_xlim(45, 45.3)
ax[1].set_ylim(-100, 1)
ax[1].grid(linestyle='--',linewidth=0.3)

plt.savefig(path_fig + 'transect_zoom_h_lat_' + str(np.round(lat[lat_index, lon_index].values,2)) + '_' + simu + '.png', dpi=300, bbox_inches='tight')

plt.close()
