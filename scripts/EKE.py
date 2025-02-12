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

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

g = xr.open_dataset(grid)
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

d = xr.open_dataset(data + simu + 'swio_avg_2017.nc')
u = d['u'][:, -1, :, :] # Vitesse surface u
v = d['v'][:, -1, :, :] # Vitesse surface v
w = d['w'][:, -1, :, :] # Vitesse surface w
d.close()

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36
u = u.where((u != fill_value), np.nan)
v = v.where((v != fill_value), np.nan)
w = w.where((w != fill_value), np.nan)

start_time = '2017-08-31'
end_time = '2017-08-31'

u_yr = u.mean(dim='time')
v_yr = v.mean(dim='time')
w_yr = w.mean(dim='time')
print('u_yr,v_yr,w_yr: ', u_yr.shape, v_yr.shape, w_yr.shape)

#u_tri = u.sel(time=slice('2017-01-01', '2017-03-31')).mean(dim='time')
#v_tri = v.sel(time=slice('2017-01-01', '2017-03-31')).mean(dim='time')
#w_tri = w.sel(time=slice('2017-01-01', '2017-03-31')).mean(dim='time')

#u_tri = u.sel(time=slice('2017-04-01', '2017-06-30')).mean(dim='time')
#v_tri = v.sel(time=slice('2017-04-01', '2017-06-30')).mean(dim='time')
#w_tri = w.sel(time=slice('2017-04-01', '2017-06-30')).mean(dim='time')

u_tri = u.sel(time=slice('2017-07-01', '2017-09-30')).mean(dim='time')
v_tri = v.sel(time=slice('2017-07-01', '2017-09-30')).mean(dim='time')
w_tri = w.sel(time=slice('2017-07-01', '2017-09-30')).mean(dim='time')

#u_tri = u.sel(time=slice('2017-10-01', '2017-12-31')).mean(dim='time')
#v_tri = v.sel(time=slice('2017-10-01', '2017-12-31')).mean(dim='time')
#w_tri = w.sel(time=slice('2017-10-01', '2017-12-31')).mean(dim='time')

u = u.sel(time=slice(start_time, end_time))
v = v.sel(time=slice(start_time, end_time))
w = w.sel(time=slice(start_time, end_time))
print('u,v,w: ', u.shape, v.shape, w.shape)

# Moyenne sur start - end
u_mean = u.mean(dim='time')
v_mean = v.mean(dim='time')
w_mean = w.mean(dim='time')

ut = u_yr - u_mean
vt = v_yr - v_mean
wt = w_yr - w_mean

print('ut,vt,wt: ', ut.shape, vt.shape, wt.shape)

# Transformation des composantes de vent (grille déformée -> grille géographique)
ut_geo = ut[:-1,:].data * np.cos(angle[:-1,:-1]) - vt[:,:-1].data * np.sin(angle[:-1,:-1])
vt_geo = ut[:-1,:].data * np.sin(angle[:-1,:-1]) + vt[:,:-1].data * np.cos(angle[:-1,:-1])
wt_geo = wt[:-1,:-1].data

print('ut,vt,wt: ', ut_geo.shape, vt_geo.shape, wt_geo.shape)
# Calcul de l'EKE
EKE = 1 / 2 * (ut_geo ** 2 + vt_geo ** 2 + wt_geo ** 2)

print('EKE: ', EKE.shape)
print('max EKE: ', np.round(float(np.max(EKE)), 3))

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_title(f"EKE SWIO {start_time}", size=9) 

cmap = cmcrameri.cm.lapaz
a = 1e-2
b = 1
c = 10
levels = np.logspace(np.log10(a), np.log10(b), c * 2 - 1) 
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon[:-1,:-1], lat[:-1,:-1], EKE, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
# Add contours / mask
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.3)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 8, 'color': 'k'}
gl.ylabel_style = {'size': 8, 'color': 'k'}

### colorbar
cb = fig.colorbar(pcm, ax=ax, label='EKE [$m^2.s^{-2}$]', norm=mpl.colors.LogNorm(vmin=a, vmax=b))
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])
cb.set_ticks(np.logspace(np.log10(a), np.log10(b), c))
cb.ax.set_yticklabels(np.round(np.logspace(np.log10(a), np.log10(b), c), 2), fontsize=8)
cb.ax.yaxis.label.set_font_properties(mpl.font_manager.FontProperties(size=8))

plt.show()

