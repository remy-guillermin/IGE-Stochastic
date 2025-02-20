import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle

# Load data (same as your existing code)
data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/'

g = xr.open_dataset(grid)
h = g['h'][:, :] # Bathymetry
lon = g['lon_rho'][:, :]  # Longitude
lat = g['lat_rho'][:, :]  # Latitude
angle = g['angle'][:, :]  # Deformation
msk = g['mask_rho'][:, :]  # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

### open netcdf file an variables
ds = xr.open_dataset(data + simu + 'swio_avg_2017.nc') # , engine='h5netcdf')
temp = ds['temp'][:, :, :, :]
s_rho = ds['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = ds['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = ds['hc'].values
ds.close()

msk_inv = np.where(msk==0, msk, np.nan)

def calc_depth(s, Cs, hc, h):
    N = len(s_rho)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h ## (hc * s[k] + h * Cs[k])
    return depth

depth_sigma = calc_depth(s_rho, Cs_rho, hc, h)

lat_index = 249
lon_index = 195

# along longitude axis
fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})

ax[0].plot(lat[:, lon_index].values, -h[:, lon_index].values, marker='o', linestyle='-', color='k', markersize=1)
ax[0].fill_between(lat[:, lon_index], -h[:, lon_index].values, y2=min(-h[:, lon_index]), color='lightgrey')

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[0].plot(lat[:, lon_index], depth_sigma[k, :, lon_index], color='grey', linestyle='-', linewidth=0.5)
ax[0].set_xlim(-16.5, 3)
ax[0].set_ylim(np.min(-h[:, lon_index].values), 5)
#ax[0].set_xlabel('Latitudes along ' + str(np.round(lon[lat_index, lon_index].values,2)) + '°E Longitude')
#ax[0].set_ylabel('Depth h (m)')
ax[0].set_xlabel('Latitudes le long de la ' + str(np.round(lon[lat_index, lon_index].values,2)) + '°E Longitude')
ax[0].set_ylabel('Profondeur h (m)')
ax[0].grid(linestyle='--',linewidth=0.3)

### subplot
ax[1].plot(lat[:, lon_index].values, -h[:, lon_index].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[1].fill_between(lat[:, lon_index], -h[:, lon_index].values, y2=min(-h[:, lon_index]), color='lightgrey')

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[1].plot(lat[:, lon_index], depth_sigma[k, :, lon_index], color='grey', linestyle='-', linewidth=0.5)
ax[1].set_xlim(-13.5, -12.5)
ax[1].set_ylim(-400, 1)
ax[1].grid(linestyle='--',linewidth=0.3)

plt.savefig(figures + 'transect_zoom_h_lon_' + str(np.round(lon[lat_index, lon_index].values,2)) + '_' + simu[:-1] + '.png', dpi=300, bbox_inches='tight')
plt.show()

# along latitude axis
fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})
ax[0].plot(lon[lat_index,:].values, -h[lat_index,:].values, marker='o', linestyle='-', color='k', markersize=1)
ax[0].fill_between(lon[lat_index,:], -h[lat_index,:].values, y2=min(-h[lat_index, :]), color='lightgrey')
# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[0].plot(lon[lat_index, :], depth_sigma[k, lat_index, :], color='grey', linestyle='-', linewidth=0.5)
ax[0].set_xlim(40, 49)
ax[0].set_ylim(-4000, 0)
#ax[0].set_xlabel('Longitudes along ' + str(np.round(lat[lat_index, lon_index].values,2)) + '°S Latitude')
#ax[0].set_ylabel('Depth (m)')
ax[0].set_xlabel('Longitudes le long de ' + str(np.round(lat[lat_index, lon_index].values,2)) + '°S Latitude')
ax[0].set_ylabel('Profondeur h (m)')
ax[0].grid(linestyle='--',linewidth=0.3)

### subplot
ax[1].plot(lon[lat_index,:].values, -h[lat_index,:].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax[1].fill_between(lon[lat_index,:], -h[lat_index,:].values, y2=min(-h[lat_index, :]), color='lightgrey')
# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax[1].plot(lon[lat_index, :], depth_sigma[k, lat_index, :], color='grey', linestyle='-', linewidth=0.5)
ax[1].set_xlim(46, 47)
ax[1].set_ylim(-500, 1)
ax[1].grid(linestyle='--',linewidth=0.3)

plt.savefig(figures + 'transect_zoom_h_lat_' + str(np.round(lat[lat_index, lon_index].values,2)) + '_' + simu[:-1] + '.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()