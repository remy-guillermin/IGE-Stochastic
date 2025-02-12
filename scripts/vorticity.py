import numpy as np
from numpy.fft import fft2, fftshift
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import sys
import cartopy.crs as ccrs
from metpy.units import units

R = 6371e3

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

g = xr.open_dataset(grid)
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
pm = g['pm'][:-1,:-1] # Curvilinear coordinate metric in XI [m^-1]
pn = g['pn'][:-1,:-1] # Curvilinear coordinate metric in ETA [m^-1]
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

d = xr.open_dataset(data + simu + 'swio_avg_2017.nc')
u = d['u'][:, -1, :, :] # Vitesse surface u
v = d['v'][:, -1, :, :] # Vitesse surface v
d.close()

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
u = u.where((u != fill_value), np.nan)
v = v.where((v != fill_value), np.nan)

start_time = '2017-07-11'
end_time = '2017-07-11'
u = u.sel(time=slice(start_time, end_time))
v = v.sel(time=slice(start_time, end_time))
print('u,v: ', u.shape, v.shape)

# Moyenne sur start - end
u_mean = u.mean(dim='time')
v_mean = v.mean(dim='time')

# Transformation des composantes de vent (grille déformée -> grille géographique)
u_geo = u_mean[:-1,:].data * np.cos(angle[:-1,:-1]) - v_mean[:,:-1].data * np.sin(angle[:-1,:-1])
v_geo = u_mean[:-1,:].data * np.sin(angle[:-1,:-1]) + v_mean[:,:-1].data * np.cos(angle[:-1,:-1])

print('u,v: ', u_geo.shape, v_geo.shape)

velocity = np.sqrt(u_geo**2 + v_geo**2)

# Calcul des dérivées
dv_dlon = np.gradient(v_geo, axis=1) * pm
du_dlat = np.gradient(u_geo, axis=0) * pn

# Calcul du rotationnel en coordonnées géographiques
vorticity = (dv_dlon - du_dlat)

helicity = velocity * vorticity

print('velocity: ', np.round(np.nanmin(velocity.values),3), np.round(np.nanmax(velocity.values),3))
print('vorticity: ', np.round(np.nanmin(vorticity.values),5), np.round(np.nanmax(vorticity.values),5))
print('helicity: ', np.round(np.nanmin(helicity.values),5), np.round(np.nanmax(helicity.values),5))

fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Define common gridline styles
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

# --- Velocity Plot ---
ax = axes[0]
ax.set_title(f"Velocity SWIO {start_time}", size=9) 

cmap = cmcrameri.cm.oslo
a, b, c = 0, 2.5, 10
levels = np.linspace(a, b, c * 2 - 1) 
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon[:-1,:-1], lat[:-1,:-1], velocity, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}

cb = fig.colorbar(pcm, ax=ax, label='Velocity [$m.s^{-1}$]', orientation='vertical')
cb.set_ticks(np.linspace(a, b, c))
cb.ax.set_yticklabels(np.round(np.linspace(a, b, c),2), fontsize=8)

# --- Vorticity Plot ---
ax = axes[1]
ax.set_title(f"Vorticity SWIO {start_time}", size=9) 

cmap = cmcrameri.cm.vik
a, b, c = -0.15, 0.15, 10
levels = np.linspace(a, b, c * 2 - 1) 
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon[:-1,:-1], lat[:-1,:-1], vorticity * 3600, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl.top_labels = False
gl.right_labels = False

cb = fig.colorbar(pcm, ax=ax, label='Vorticity [$h^{-1}$]', orientation='vertical')
cb.set_ticks(np.linspace(a, b, c))
cb.ax.set_yticklabels(np.round(np.linspace(a, b, c),2), fontsize=8)

# --- Helicity Plot ---
ax = axes[2]
ax.set_title(f"Helicity SWIO {start_time}", size=9) 

cmap = cmcrameri.cm.vik
a, b, c = -0.5, 0.5, 11
levels = np.linspace(a, b, c * 2 - 1) 
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon[:-1,:-1], lat[:-1,:-1], helicity * 3600 ** 2 / 1000, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl.top_labels = False
gl.right_labels = False

cb = fig.colorbar(pcm, ax=ax, label='Helicity [$km.h^{-2}$]', orientation='vertical')
cb.set_ticks(np.linspace(a, b, c))
cb.ax.set_yticklabels(np.round(np.linspace(a, b, c),2), fontsize=8)

plt.tight_layout()
plt.show()


