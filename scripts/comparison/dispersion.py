# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmocean
import cartopy.crs as ccrs
import glob
import os
from functools import partial



figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio.nc'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

SWIO = (25, 69, -36, 7)
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}


# Replace "NaN_MERGED" with the folder you want to use
# Options are "NaN_MERGED", "MLD" or "WINDSTR"

ini_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_ini', 'mean*.nc'))[0]
print(f'ini mean file found: {ini_mean_file}')
ini_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_ini', 'std*.nc'))[0]
print(f'ini std file found: {ini_std_file}')

str_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_str', 'mean*.nc'))[0]
print(f'str mean file found: {str_mean_file}')
str_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_str', 'std*.nc'))[0]
print(f'str std file found: {str_std_file}')

gls_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_gls', 'mean*.nc'))[0]
print(f'gls mean file found: {gls_mean_file}')
gls_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_gls', 'std*.nc'))[0]
print(f'gls std file found: {gls_std_file}')


ini_mean = xr.open_dataset(ini_mean_file).sortby('time')
ini_std = xr.open_dataset(ini_std_file).sortby('time')

str_mean = xr.open_dataset(str_mean_file).sortby('time')
str_std = xr.open_dataset(str_std_file).sortby('time')

gls_mean = xr.open_dataset(gls_mean_file).sortby('time')
gls_std = xr.open_dataset(gls_std_file).sortby('time')

g = xr.open_dataset(grid)
lon = g.lon_rho
lat = g.lat_rho
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')



date = str(ini_mean.time.isel(time=-100).data.astype('datetime64[D]'))
print(f"Plotting {date}.")

fig, axs = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()}, sharey=True)
fig.suptitle(f"Ensemble salinity standard deviation {date}")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black')
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}
    

    
ax = axs[0]
ax.set_title('Initial condition')

pcm = ax.pcolormesh(lon, lat, ini_std.salt.sel(time=date).isel(time=0), cmap='magma', vmin=0, vmax=0.4)    
    
ax = axs[1]
ax.set_title('Drag coefficient')   

pcm = ax.pcolormesh(lon, lat, str_std.salt.sel(time=date).isel(time=0), cmap='magma', vmin=0, vmax=0.4)   
    
ax = axs[2]
ax.set_title('Vertical mixing') 

pcm = ax.pcolormesh(lon, lat, gls_std.salt.sel(time=date).isel(time=0), cmap='magma', vmin=0, vmax=0.4)  

cb = plt.colorbar(pcm, ax=ax, orientation='vertical', pad=0.05, label='Standard deviation [psu]')   
    
fig.tight_layout()
fig.savefig(f'{figures}/comp_std_{date}.png', dpi=300, transparent=True)
plt.show()
