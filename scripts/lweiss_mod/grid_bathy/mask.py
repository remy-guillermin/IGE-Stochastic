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
data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/RUN_CROCO/'
grid = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
simu = 'run_swio_stogen_1mth_diff_100p_10jo1ud12_0001'

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

msk_inv = np.where(msk==0, msk, np.nan)

#sys.exit()

### open netcdf file an variables
ds = xr.open_dataset(data + simu + '/001swiose_his.nc') # , engine='h5netcdf')
temp = ds['temp'][:, :, :, :]
#h = ds['h'][:, :]  # bathymetry at RHO-points
#lon = ds['lon_rho'][:, :]
#lat = ds['lat_rho'][:, :]
#msk = ds['mask_rho'][:, :]
s_rho = ds['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = ds['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = ds['hc'].values # S-coordinate parameter, critical depth
ds.close()

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
# Coordonnées des points à zoomer
point = (249, 195)
# point = (249, 195)
# Définition des indices de grille pour le zoom
S = 40
zoom = (slice(point[0]-int(S/2), point[0]+int(S/2)), slice(point[1]-int(S/2), point[1]+int(S/2)))

#######################################################################################
### plot
#######################################################################################
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

cs = ax.contour(lon, lat, h, levels=np.arange(0, 6000, 1000), colors='grey',
                linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs, fmt='%d', inline=True, fontsize=4)
cs2 = ax.contour(lon, lat, h, levels=np.arange(-300, 600, 400), colors='grey', 
                linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs2, fmt='%d', inline=True, fontsize=4)
cs3 = ax.contour(lon, lat, h, levels=np.arange(0, 100, 50), colors='red', linewidths=0.3)
plt.clabel(cs3, fmt='%d', inline=True, fontsize=5)
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 6, 'color': 'k'}
gl.ylabel_style = {'size': 6, 'color': 'k'}
#ax.coastlines(resolution='10m')
#ax.set_title('Sea Surface Temperature')

plt.show()