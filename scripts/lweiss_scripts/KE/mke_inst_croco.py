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
from metpy.units import units
import metpy.calc as mpcalc
##

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_ref2_2017'

### open netcdf file an variables
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
msk_inv = np.where(msk==0, msk, np.nan)
dx = 1 / g['pm'] * units.meter
dy = 1 / g['pn'] * units.meter
g.close()

dx = units.Quantity(dx.data, "m")
dy = units.Quantity(dy.data, "m")

ds = xr.open_dataset(data + simu + '/swio_avg.nc') # , engine='h5netcdf')
u = ds['u'][:, -1, :, :]  # Sélectionne le premier niveau s_rho
v = ds['v'][:, -1, :, :]
# ds.close()

uu = units.Quantity(u[:,:-1,:].data, "m/s")
vv = units.Quantity(v[:,:,:-1].data, "m/s")

mke = 0.5 * (uu ** 2 + vv ** 2)

#vorticities = []
#for t in range(u.shape[0]):
#    vorticity = mpcalc.vorticity(uu[t,:,:], vv[t,:,:], dx=dx[:-1,:-2], dy=dy[:-2,:-1])
#    vorticities.append(vorticity)

print('U, V, mke, lon, lat, msk \n', u.shape, v.shape, mke.shape, lon.shape, lat.shape, msk.shape)

### plot
for t in range(0, mke.shape[0]):
    # Convertir numpy.datetime64 en chaîne de caractères
    time_val = ds.time[t].values
    time_str = np.datetime_as_string(time_val, unit='s')
    # Formater la chaîne de caractères pour l'utiliser dans le titre et le nom du fichier
    t_title = time_str.replace("T", " ")
    time_str_formatted = time_str.replace("T", "-").replace(":", "-")

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    cmap = cmcrameri.cm.roma_r
    # a = 0
    # b = 0.5
    # c = 11
    # levels = np.linspace(a, b, c*2-1)  ### amplitude de la colorbar à ajuster ici
    norm = mpl.colors.LogNorm(vmin=1e-2, vmax=1e0)
    # norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    pcm = ax.pcolormesh(lon, lat, mke[t], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
    ax.contourf(lon, lat, msk_inv, colors='lightgray')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'k'}
    gl.ylabel_style = {'size': 8, 'color': 'k'}
    #ax.coastlines(resolution='10m')
    ax.set_title(f"Daily MKE {t_title[:10]}", size=10)

    ### colorbar
    cb = fig.colorbar(pcm, ax=ax, label='Mean Kinetic Energy $[m^{2} s^{-2}]$')
    # colorbaryticks = np.linspace(a, b, c)  # levels[::2]
    posax = ax.get_position()
    poscb = cb.ax.get_position()
    cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
    # cb.set_ticks(colorbaryticks)
    # cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),1), fontsize=8)
    text = cb.ax.yaxis.label
    font = mpl.font_manager.FontProperties(size=8)
    text.set_font_properties(font)

    ### SAVE
    plt.savefig(path_fig + simu + f'/mke/mke_{simu[:13]}_{time_str_formatted[:10]}.png', dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()
