import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import sys
import cartopy.crs as ccrs
from scipy.interpolate import interp1d
from metpy.units import units

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

g = xr.open_dataset(grid)
h = g['h'][:, :]
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

d = xr.open_dataset(data + simu + 'swio_avg.nc')
s_rho = d['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = d['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = d['hc'].values # S-coordinate parameter, critical depth
d.close()

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

fig = plt.figure(figsize=(6, 6))
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

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
cmap = cmocean.cm.deep_r
a = 0
b = 5000
c = int((b-a)*2/1000+1)
levels = np.linspace(a, b, c*2-1) 
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

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 6, 'color': 'k'}
gl.ylabel_style = {'size': 6, 'color': 'k'}
step = 10
for xi in range(0, lon.shape[1], step):  # Loop through the xi (longitude) points
    ax.plot(lon[:, xi], lat[:, xi], color='gray', linewidth=0.3, alpha=0.6, transform=ccrs.PlateCarree())  # Vertical lines
for eta in range(0, lat.shape[0], step):  # Loop through the eta (latitude) points
    ax.plot(lon[eta, :], lat[eta, :], color='gray', linewidth=0.3, alpha=0.6, transform=ccrs.PlateCarree())  # Horizontal lines

cb = fig.colorbar(pcm, ax=ax, label='Bathymetry (m)')
colorbaryticks = np.linspace(a, b, c) 
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=6)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=6)
text.set_font_properties(font)

plt.show()