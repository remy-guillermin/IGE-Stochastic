import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import sys
import cartopy.crs as ccrs
from metpy.units import units

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

simu = 'run_swio_stogen_1mth_diff_100p_10jo1ud12'

g = xr.open_dataset(data + simu + '_0001/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk == 0, msk, np.nan)
# dx = 1 / g['pm'] * units.meter
# dy = 1 / g['pn'] * units.meter
angle = g['angle'][:, :]
#sys.exit()
g.close()

ds = xr.open_dataset(data + simu + '_0001/001swiose_his.nc')  # , engine='h5netcdf')
u = ds['sustr'][:, :, :]
v = ds['svstr'][:, :, :] 
ds.close()

start_time = '2017-01-31' # '2019-01-01T12:00:00'
end_time = '2017-01-31'
u = u.sel(time=slice(start_time, end_time))
v = v.sel(time=slice(start_time, end_time))
# dx = units.Quantity(dx.data, "m")
# dy = units.Quantity(dy.data, "m")
print('u,v: ', u.shape, v.shape)
# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
u = u.where((u != fill_value), np.nan)
v = v.where((v != fill_value), np.nan)
# dx = dx.where(msk == 1, np.nan)
# dy = dy.where(msk == 1, np.nan)
print('u,v: ', u.shape, v.shape)
# Moyenne temporelle
u_mean = u.mean(dim='time')
v_mean = v.mean(dim='time')
# Transformation des composantes de vent (grille déformée -> grille géographique)
u_geo = u_mean[:-1,:].data * np.cos(angle[:-1,:-1]) - v_mean[:,:-1].data * np.sin(angle[:-1,:-1])
v_geo = u_mean[:-1,:].data * np.sin(angle[:-1,:-1]) + v_mean[:,:-1].data * np.cos(angle[:-1,:-1])
print('u,v: ', u_geo.shape, v_geo.shape)
# Calcul de l'intensité du stress du vent
wind_stress = np.sqrt(u_geo**2 + v_geo**2)
print('wind stress: ', wind_stress.shape)

# sys.exit()

# Tracer la figure représentant la moyenne temporelle
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_title(f"Wind Stress SWIO {start_time}", size=9) 

cmap = cmcrameri.cm.batlow #roma_r #batlowW #  plt.cm.seismic
a = 0
b = 0.2
c = 5
levels = np.linspace(a, b, c * 4 - 3)  ### amplitude de la colorbar à ajuster ici
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon[:-1,:-1], lat[:-1,:-1], wind_stress, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
# Ajouter les vecteurs de stress du vent (u_geo, v_geo)
step = 10
quiver = ax.quiver(lon[:-1:step,:-1:step].data, lat[:-1:step,:-1:step].data, u_geo[::step, ::step], v_geo[::step, ::step], transform=ccrs.PlateCarree(), width=0.05, units='xy', scale=0.09, headlength=4, headwidth=4, color='black')
# Add contours / mask
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.3)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='wind stress (N/m²)')
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])
cb.set_ticks(np.linspace(a, b, c))
cb.ax.set_yticklabels(np.round(np.linspace(a, b, c),2), fontsize=8)
cb.ax.yaxis.label.set_font_properties(mpl.font_manager.FontProperties(size=8))

plt.show()
