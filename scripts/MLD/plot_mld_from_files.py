import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import glob
%matplotlib inline



figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio.nc'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

lat_index = [430, 280]

colors = ['black', 'royalblue', 'crimson']



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
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')



fig, ax = plt.subplots(1, 1, figsize=(8,3))

time_id = 150

date = str(gls_mean.isel(time=time_id).time.values.astype('datetime64[D]'))

fig.suptitle(rf'Zonal slice MLD - {date}')

land_color = ccrs.cartopy.feature.COLORS['land']

ax.fill_between(lon[lat_index[0], :], -h[lat_index[0], :].values, y2=min(-h[lat_index[0], :]), color=land_color)
ax.plot(lon[lat_index[0], :], -h[lat_index[0], :].values, color='black', linewidth=0.5)

for ds in ini_ds:
    data = ds.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
    ax.plot(lon.isel(eta_rho=lat_index[0]), data, color=colors[0], linewidth=0.3, alpha=0.3)

mean_data = ini_mean.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
ax.plot(mean_data.lon_rho, mean_data, color=colors[0], label='INI')

for ds in str_ds:
    data = ds.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
    ax.plot(lon.isel(eta_rho=lat_index[0]), data, color=colors[1], linewidth=0.3, alpha=0.3)

mean_data = str_mean.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
ax.plot(mean_data.lon_rho, mean_data, color=colors[1], label='STR')

for ds in gls_ds:
    data = ds.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
    ax.plot(lon.isel(eta_rho=lat_index[0]), data, color=colors[2], linewidth=0.3, alpha=0.3)

mean_data = gls_mean.mld.sel(time=date).isel(time=0, eta_rho=lat_index[0])
ax.plot(mean_data.lon_rho, mean_data, color=colors[2], label='GLS')

ticks = ax.get_yticks()
ax.set_yticks(ticks)
ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
ax.grid(alpha=0.5)

ticks = ax.get_xticks()
ax.set_xticks(ticks)
ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])

ax.set_xlim((40,65))
ax.set_ylim(-110, 0)

ax.legend(loc='lower right')

fig.tight_layout()
plt.show()






