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
%matplotlib inline



figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio.nc'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

lat_index = [430, 280]

colors = ['black', 'royalblue', 'crimson']

SWIO = (25, 69, -36, 7)
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}



ini_files = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', '0*.nc'))
ini_files.sort()
print(f'{len(ini_files)} ini files found')

ini_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'mean*.nc'))[0]
print(f'ini mean file found: {ini_mean_file}')
ini_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'std*.nc'))[0]
print(f'ini std file found: {ini_std_file}')

str_files = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', '0*.nc'))
str_files.sort()
print(f'{len(str_files)} str files found')

str_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'mean*.nc'))[0]
print(f'str mean file found: {str_mean_file}')
str_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'std*.nc'))[0]
print(f'str std file found: {str_std_file}')


gls_files = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', '0*.nc'))
gls_files.sort()
print(f'{len(gls_files)} gls files found')

gls_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'mean*.nc'))[0]
print(f'gls mean file found: {gls_mean_file}')
gls_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'std*.nc'))[0]
print(f'gls std file found: {gls_std_file}')



ini_ds = [xr.open_dataset(f) for f in ini_files]
ini_mean = xr.open_dataset(ini_mean_file)
ini_std = xr.open_dataset(ini_std_file)

str_ds = [xr.open_dataset(f) for f in str_files]
str_mean = xr.open_dataset(str_mean_file)
str_std = xr.open_dataset(str_std_file)

gls_ds = [xr.open_dataset(f) for f in gls_files]
gls_mean = xr.open_dataset(gls_mean_file)
gls_std = xr.open_dataset(gls_std_file)

g = xr.open_dataset(grid)
lon = g.lon_rho
lat = g.lat_rho
xi = g.xi_rho
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')



fig, ax = plt.subplots(1, 1, figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"Slices")

ax.set_extent(SWIO)
ax.coastlines(resolution='50m')
ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=1)
ax.add_feature(ccrs.cartopy.feature.OCEAN, zorder=0)
ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=1)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=1)

land_color = ccrs.cartopy.feature.COLORS['land']

minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
    category='physical',
    name='minor_islands',
    scale='10m',
    facecolor=land_color,
    edgecolor='black'
)

ax.add_feature(minor_islands, zorder=1)

ax.plot(lon[lat_index[0], :], lat[lat_index[0], :], color='red', linewidth=2, zorder=0)
ax.plot(lon[lat_index[1], :], lat[lat_index[1], :], color='red', linewidth=2, zorder=0)

ax.text(lon[lat_index[0], :].values[-1] + 0.5, lat[lat_index[0], :].values[-1] - 0.5, 'S1', fontsize=12, color='red')
ax.text(lon[lat_index[1], :].values[-1] + 0.5, lat[lat_index[1], :].values[-1] - 0.5, 'S2', fontsize=12, color='red')


gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = gl.ylabel_style = {'color': 'k'}



fig.tight_layout()
fig.savefig(os.path.join(figures, 'slice.png'), dpi=300, transparent=True)
plt.show()



time_id = 250 #100, 160, 620
fig, axs = plt.subplots(3, 2, figsize=(16,9), sharex=True, sharey=True)
# fig, axs = plt.subplots(2, 2, figsize=(16,6), sharex=True, sharey=True)

date = str(gls_mean.isel(time=time_id).time.values.astype('datetime64[D]'))

fig.suptitle(rf'Zonal slice MLD - {date}')

land_color = ccrs.cartopy.feature.COLORS['land']

axs[0,0].set_title('Slice S1')
axs[0,1].set_title('Slice S2')

for ax in axs[:,0]:        
    ax.fill_between(lon[lat_index[0], :], -h[lat_index[0], :].values, y2=min(-h[lat_index[0], :]), color=land_color)
    ax.plot(lon[lat_index[0], :], -h[lat_index[0], :].values, color='black', linewidth=0.5)
    
    ax.set_xlim((40,65))
    ax.set_ylim(-110, 0)
    
for ax in axs[:,1]:        
    ax.fill_between(lon[lat_index[1], :], -h[lat_index[1], :].values, y2=min(-h[lat_index[1], :]), color=land_color)
    ax.plot(lon[lat_index[1], :], -h[lat_index[1], :].values, color='black', linewidth=0.5)
    
    ax.set_xlim((40,65))
        
for i, ax in enumerate(axs[0,:]):
    id = lat_index[i]
    for ds in ini_ds:
        data = ds.mld.sel(time=date).isel(time=0, eta_rho=id)
        ax.plot(lon.isel(eta_rho=id), data, color=colors[0], linewidth=0.3, alpha=0.5)
        
    ax.plot(lon.isel(eta_rho=id), data, color=colors[0], linewidth=0.3, alpha=0.5, label='Member')
    
    mean_data = ini_mean.mld.sel(time=date).isel(time=0, eta_rho=id)
    rms = np.sqrt(np.nanmean((ini_std.mld.sel(time=date).isel(eta_rho=id) ** 2)))
    ax.plot(mean_data.lon_rho, mean_data, color=colors[0], label=f'Ensemble mean')
    ax.legend(
        title=f'INI - RMSE: {rms:.1f} m', 
        loc='lower right', framealpha=1.0, 
        edgecolor='black', 
        title_fontproperties={'weight':'semibold'}
    )
 
for i, ax in enumerate(axs[1,:]):
    id = lat_index[i]
    for ds in str_ds:
        data = ds.mld.sel(time=date).isel(time=0, eta_rho=id)
        ax.plot(lon.isel(eta_rho=id), data, color=colors[1], linewidth=0.3, alpha=0.5)

    ax.plot(lon.isel(eta_rho=id), data, color=colors[1], linewidth=0.3, alpha=0.5, label='Member')
    

    mean_data = str_mean.mld.sel(time=date).isel(time=0, eta_rho=id)
    rms = np.sqrt(np.nanmean((str_std.mld.sel(time=date).isel(eta_rho=id) ** 2)))
    ax.plot(mean_data.lon_rho, mean_data, color=colors[1], label=f'Ensemble mean')
    ax.legend(
        title=f'STR - RMSE: {rms:.1f} m', 
        loc='lower right', framealpha=1.0, 
        edgecolor='black',
        title_fontproperties={'weight':'semibold'}
        )
 
for i, ax in enumerate(axs[2,:]):
    id = lat_index[i]
    for ds in gls_ds:
        data = ds.mld.sel(time=date).isel(time=0, eta_rho=id)
        ax.plot(lon.isel(eta_rho=id), data, color=colors[2], linewidth=0.3, alpha=0.5)

    ax.plot(lon.isel(eta_rho=id), data, color=colors[2], linewidth=0.3, alpha=0.5, label='Member')
    

    mean_data = gls_mean.mld.sel(time=date).isel(time=0, eta_rho=id)
    rms = np.sqrt(np.nanmean((gls_std.mld.sel(time=date).isel(eta_rho=id) ** 2)))
    ax.plot(mean_data.lon_rho, mean_data, color=colors[2], label=f'Ensemble mean')
    ax.legend(
        title=f'GLS - RMSE: {rms:.1f} m', 
        loc='lower right', 
        framealpha=1.0, edgecolor='black', 
        title_fontproperties={'weight':'semibold'}
    )

for ax in axs[:,0]:
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
for ax in axs.flatten():
    ax.grid(alpha=0.5)
    
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])

fig.tight_layout()
fig.savefig(f'{figures}/slices_mld_comp_{date}.png', dpi=300, transparent=True)
plt.show()



std = (ini_std.mld.sel(time=date).isel(eta_rho=lat_index[1], time=0) ** 2)
rms_left = np.sqrt(np.nanmean(std.isel(xi_rho=slice(0,242))))
rms_right = np.sqrt(np.nanmean(std.isel(xi_rho=slice(242,-1))))
rms_left, rms_right



std = (str_std.mld.sel(time=date).isel(eta_rho=lat_index[1], time=0) ** 2)
rms_left = np.sqrt(np.nanmean(std.isel(xi_rho=slice(0,242))))
rms_right = np.sqrt(np.nanmean(std.isel(xi_rho=slice(242,-1))))
rms_left, rms_right



std = (gls_std.mld.sel(time=date).isel(eta_rho=lat_index[1], time=0) ** 2)
rms_left = np.sqrt(np.nanmean(std.isel(xi_rho=slice(0,242))))
rms_right = np.sqrt(np.nanmean(std.isel(xi_rho=slice(242,-1))))
print(rms_left, rms_right)