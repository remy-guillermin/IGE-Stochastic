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

def bvor(ssh,dx,dy,lat):
    # Calcul des dérivées de ssh
    dssh_dy = np.gradient(ssh, axis=0) / dy
    dssh_dx = np.gradient(ssh, axis=1) / dx
    # Calcul de f (paramètre de Coriolis)
    # Omega = 7.2921e-5  # s^-1, vitesse angulaire de la Terre
    # lat_rad = np.deg2rad(lat)
    # f = 2 * Omega * np.sin(lat_rad)
    # Calcul des vitesses géostrophiqueq
    g = 9.81  # m/s^2
    u = -(g / f) * dssh_dy
    v =  (g / f) * dssh_dx
    # Calcul des dérivées de u et v
    du_dy = np.gradient(u, axis=0) / dy
    dv_dx = np.gradient(v, axis=1) / dx
    # Calcul de la vorticité relative
    vorticity = (dv_dx - du_dy)
    return vorticity

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_0003_2017'

### open grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
f = g['f'][:,:] # Coriolis parameter
dx = 1 / g['pm'] 
dy = 1 / g['pn']
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk == 0, msk, np.nan)
#sys.exit()
g.close()

ds = xr.open_dataset(data + simu + '/swiose_avg.nc')  # , engine='h5netcdf')
zeta = ds['zeta'][:, :, :]  # Sélectionne le premier niveau s_rho
# lon = ds['lon_rho'][:, :]
# lat = ds['lat_rho'][:, :]
# msk = ds['mask_rho'][:, :]

# Print the shape for debugging purposes
print(zeta.shape, lon.shape, lat.shape, msk.shape)

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
zeta = zeta.where(zeta != fill_value, np.nan)
dx = dx.where(msk == 1, np.nan)
dy = dy.where(msk == 1, np.nan)

# Définir la période pour calculer la moyenne temporelle
start_time = '2017-08-01' # '2019-01-01T12:00:00'
end_time = '2017-08-21'

# Sélectionner la période
zeta_period = zeta.sel(time=slice(start_time, end_time))

vorticities = []
for t in range(zeta_period.shape[0]):
    vorticity = bvor(zeta_period[t,:,:], dx, dy, lat)
    vorticities.append(vorticity)

# Calculer la moyenne temporelle sur la période sélectionnée
vorticity_mean = np.mean(np.array(vorticities), axis=0) # vorticity.mean(dim='time', skipna=True)

# Tracer la figure représentant la moyenne temporelle
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

cmap = plt.cm.seismic
a = -0.5
b = 0.5
c = 11 # int((b - a) + 1)
levels = np.linspace(a, b, c * 2 - 1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, vorticity_mean * 10000, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}
# ax.coastlines(resolution='10m')
ax.set_title(f"Surface vorticity SWIO from {start_time} to {end_time}", size=9)

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='vorticity $[10^4 s^{-1}]$')
colorbaryticks = np.linspace(a, b, c)  # levels[::2]
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),1), fontsize=8)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=8)
text.set_font_properties(font)

### SAVE
plt.savefig(path_fig + simu + f'/Vorticity/vorticity_mean_{simu}_{start_time}_{end_time}.png', dpi=300, bbox_inches='tight')
plt.close()

ds.close()

