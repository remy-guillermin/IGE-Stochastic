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
import metpy.calc as mpcalc
from metpy.units import units

data = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

simu = 'run_swio2_deter_2017_2023'

### open netcdf file and variables
### open grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk == 0, msk, np.nan)
dx = 1 / g['pm'] * units.meter
dy = 1 / g['pn'] * units.meter
#sys.exit()
g.close()

ds = xr.open_dataset(data + simu + '/swio_avg.nc')  # , engine='h5netcdf')
u = ds['u'][:, -1, :, :]  # Sélectionne le premier niveau s_rho
v = ds['v'][:, -1, :, :] 
# lon = ds['lon_rho'][:, :]
# lat = ds['lat_rho'][:, :]
# msk = ds['mask_rho'][:, :]

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
u = u.where((u != fill_value), np.nan)
v = v.where((v != fill_value), np.nan)
dx = dx.where(msk == 1, np.nan)
dy = dy.where(msk == 1, np.nan)

# Définir la période pour calculer la moyenne temporelle
start_time = '2017-10-01' # '2019-01-01T12:00:00'
end_time = '2017-10-02'

# Sélectionner la période
u_period = u.sel(time=slice(start_time, end_time))
v_period = v.sel(time=slice(start_time, end_time))

dx = units.Quantity(dx.data, "m")
dy = units.Quantity(dy.data, "m")

# dx = dx.expand_dims(dim={'time': u_period.shape[0]}, axis=0) 
# dy = dy.expand_dims(dim={'time': u_period.shape[0]}, axis=0) 

# lon = lon.expand_dims(dim={'time': u_period.shape[0]}, axis=0)
# lat = lat.expand_dims(dim={'time': u_period.shape[0]}, axis=0)

# u_period = u_period.rename({'xi_u':'xi_rho'})
# v_period = v_period.rename({'eta_v':'eta_rho'})

uu = units.Quantity(u_period[:,:-1,:].data, "m/s")
vv = units.Quantity(v_period[:,:,:-1].data, "m/s") 

vorticities = []
for t in range(u_period.shape[0]):
    vorticity = mpcalc.vorticity(uu[t,:,:], vv[t,:,:], dx=dx[:-1,:-2], dy=dy[:-2,:-1])
    vorticities.append(vorticity)

# Print the shape for debugging purposes
print(vorticity.shape, lon.shape, lat.shape, msk.shape)

# Calculer la moyenne temporelle sur la période sélectionnée
vort_mean = np.mean(np.array(vorticities), axis=0) # vorticity.mean(dim='time', skipna=True)

print(vort_mean.shape)

# Tracer la figure représentant la moyenne temporelle
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

cmap = plt.cm.seismic
a = -0.5
b = 0.5
c = 11
levels = np.linspace(a, b, c * 2 - 1)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, vort_mean * 10000, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
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
cb = fig.colorbar(pcm, ax=ax, label='vorticity $[10^{4} s^{-1}]$')
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
plt.show()

ds.close()

