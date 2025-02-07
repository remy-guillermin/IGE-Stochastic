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
# Coordonnées des points à zoomer
point = (249, 195)
# point = (249, 195)
# Définition des indices de grille pour le zoom
S = 40
zoom = (slice(point[0]-int(S/2), point[0]+int(S/2)), slice(point[1]-int(S/2), point[1]+int(S/2)))

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

#######################################################################################
#### h transect
#######################################################################################
lat_index = 243 # 249
lon_index = 162 # 195

# along longitude axis
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(lat[:, lon_index].values, -h[:, lon_index].values, marker='o', linestyle='-', color='k',
        markersize=1)
ax.fill_between(lat[:, lon_index], -h[:, lon_index].values, y2=min(-h[:, lon_index]), color='lightgrey')
# ax.contourf(lat[:, lon_index].values, -h[:, lon_index].values, color='lightgrey', alpha=0.5)

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax.plot(lat[:, lon_index], depth_sigma[k, :, lon_index], color='grey', linestyle='-', linewidth=0.5)
# ax.scatter(lat[lat_index,lon_index].values, -h[lat_index,lon_index].values, color='red', marker='o', s=50)
ax.set_xlim(-16.5, 3)
ax.set_ylim(np.min(-h[:, lon_index].values), 5)
ax.set_xlabel('Latitudes along ' + str(np.round(lon[lat_index, lon_index].values,2)) + '°E Longitude')
ax.set_ylabel('Depth h (m)')
plt.grid(linestyle='--',linewidth=0.3)
plt.show()

#######################################################################################
# along latitude axis
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(lon[lat_index,:].values, -h[lat_index,:].values, marker='o', linestyle='-', color='k',
        markersize=1)

# Ajouter des courbes de niveau
for k in range(len(s_rho)):
    ax.plot(lon[lat_index, :], depth_sigma[k, lat_index, :], color='grey', linestyle='--', linewidth=0.5)

# ax.scatter(lon[lat_index,lon_index].values, -h[lat_index,lon_index].values, color='red', marker='o', s=50)
ax.set_xlim(39, 49)
ax.set_ylim(-4000, 0)
ax.set_xlabel('Longitudes along ' + str(np.round(lat[lat_index, lon_index].values,2)) + '°S Latitude')
ax.set_ylabel('Depth h (m)')
plt.grid()
plt.show()
