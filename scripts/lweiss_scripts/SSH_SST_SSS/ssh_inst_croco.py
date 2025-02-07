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
##

data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_ref2_2017'

# boxes = [(35,40,-23,-19)] # MoC

### open netcdf file an variables
### open grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
#sys.exit()
g.close()

ds = xr.open_dataset(data + simu + '/swio_avg.nc') # , engine='h5netcdf')
zeta = ds['zeta'][:, :, :]  # Sélectionne le premier niveau s_rho
# lon = ds['lon_rho'][:,:]
# lat = ds['lat_rho'][:,:]
# msk = ds['mask_rho'][:,:]
# ds.close()

msk_inv = np.where(msk==0, msk, np.nan)

print(zeta.shape, lon.shape, lat.shape, msk.shape)

### plot
for t in range(0,zeta.shape[0]):
    # Convertir numpy.datetime64 en chaîne de caractères
    time_val = ds.time[t].values
    time_str = np.datetime_as_string(time_val, unit='s')
    # Formater la chaîne de caractères pour l'utiliser dans le titre et le nom du fichier
    t_title = time_str.replace("T", " ")
    time_str_formatted = time_str.replace("T", "-").replace(":", "-")

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    cmap = cmcrameri.cm.roma_r
    a = 0
    b = 1.2
    c = 13 # int((b-a)+1)
    levels = np.linspace(a, b, c)  ### amplitude de la colorbar à ajuster ici
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    pcm = ax.pcolormesh(lon, lat, zeta[t], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
    ax.contourf(lon, lat, msk_inv, colors='lightgray')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'k'}
    gl.ylabel_style = {'size': 8, 'color': 'k'}
    #ax.coastlines(resolution='10m')
    ax.set_title(f"Daily mean SSH {t_title[:10]}", size=10)

    ### colorbar
    cb = fig.colorbar(pcm, ax=ax, label='Sea Surface Height [m]')
    colorbaryticks = np.linspace(a, b, c)  # levels[::2]
    posax = ax.get_position()
    poscb = cb.ax.get_position()
    cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])  # [poscb.x0, posax.y0, poscb.width, posax.height]
    cb.set_ticks(colorbaryticks)
    cb.ax.set_yticklabels(np.round(colorbaryticks.astype(float),2), fontsize=8)
    text = cb.ax.yaxis.label
    font = mpl.font_manager.FontProperties(size=8)
    text.set_font_properties(font)

    ### SAVE
    plt.savefig(path_fig + simu + f'/ssh/ssh_{simu[:13]}_{time_str_formatted[:10]}.png', dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()
# sys.exit()
