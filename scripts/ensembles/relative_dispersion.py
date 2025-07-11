# LIBRAIRIES
import numpy as np
import xarray as xr
import pandas as pd
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

lat_index = [430, 280]

ensembles = ['GLS', 'STR', 'INI']
ens_colors = ['crimson', 'royalblue', 'black']
ens_linestyles = ['solid', 'dashed', 'dotted']

box_colors = ['sandybrown', 'slateblue', 'lightseagreen', 'red']
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']

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


ini_mean_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'mean*.nc'))[0]
print(f'ini mld mean file found: {ini_mean_mld_file}')
ini_std_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'std*.nc'))[0]
print(f'ini mld std file found: {ini_std_mld_file}')

str_mean_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'mean*.nc'))[0]
print(f'str mld mean file found: {str_mean_mld_file}')
str_std_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'std*.nc'))[0]
print(f'str mld std file found: {str_std_mld_file}')

gls_mean_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'mean*.nc'))[0]
print(f'gls mld mean file found: {gls_mean_mld_file}')
gls_std_mld_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'std*.nc'))[0]
print(f'gls mld std file found: {gls_std_mld_file}')


ini_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_ini', 'mean*.nc'))[0]
print(f'ini str mean file found: {ini_mean_str_file}')
ini_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_ini', 'std*.nc'))[0]
print(f'ini str std file found: {ini_std_str_file}')

str_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_str', 'mean*.nc'))[0]
print(f'str str mean file found: {str_mean_str_file}')
str_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_str', 'std*.nc'))[0]
print(f'str str std file found: {str_std_str_file}')

gls_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'mean*.nc'))[0]
print(f'gls str mean file found: {gls_mean_str_file}')
gls_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'std*.nc'))[0]
print(f'gls str std file found: {gls_std_str_file}')

ini_mean = xr.open_dataset(ini_mean_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
ini_std = xr.open_dataset(ini_std_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

ini_mean = ini_mean.drop_duplicates(dim='time', keep='last')
ini_std = ini_std.drop_duplicates(dim='time', keep='last')

str_mean = xr.open_dataset(str_mean_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
str_std = xr.open_dataset(str_std_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

str_mean = str_mean.drop_duplicates(dim='time', keep='last')
str_std = str_std.drop_duplicates(dim='time', keep='last')

gls_mean = xr.open_dataset(gls_mean_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
gls_std = xr.open_dataset(gls_std_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

gls_mean = gls_mean.drop_duplicates(dim='time', keep='last')
gls_std = gls_std.drop_duplicates(dim='time', keep='last')

ini_mld_mean = xr.open_dataset(ini_mean_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
ini_mld_std = xr.open_dataset(ini_std_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

ini_mld_mean = ini_mld_mean.drop_duplicates(dim='time', keep='last')
ini_mld_std = ini_mld_std.drop_duplicates(dim='time', keep='last')

str_mld_mean = xr.open_dataset(str_mean_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
str_mld_std = xr.open_dataset(str_std_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

str_mld_mean = str_mld_mean.drop_duplicates(dim='time', keep='last')
str_mld_std = str_mld_std.drop_duplicates(dim='time', keep='last')

gls_mld_mean = xr.open_dataset(gls_mean_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
gls_mld_std = xr.open_dataset(gls_std_mld_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

gls_mld_mean = gls_mld_mean.drop_duplicates(dim='time', keep='last')
gls_mld_std = gls_mld_std.drop_duplicates(dim='time', keep='last')

ini_str_mean = xr.open_dataset(ini_mean_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
ini_str_std = xr.open_dataset(ini_std_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

ini_str_mean = ini_str_mean.drop_duplicates(dim='time', keep='last')
ini_str_std = ini_str_std.drop_duplicates(dim='time', keep='last')

str_str_mean = xr.open_dataset(str_mean_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
str_str_std = xr.open_dataset(str_std_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

str_str_mean = str_str_mean.drop_duplicates(dim='time', keep='last')
str_str_std = str_str_std.drop_duplicates(dim='time', keep='last')

gls_str_mean = xr.open_dataset(gls_mean_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))
gls_str_std = xr.open_dataset(gls_std_str_file).sortby('time').isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))

gls_str_mean = gls_str_mean.drop_duplicates(dim='time', keep='last')
gls_str_std = gls_str_std.drop_duplicates(dim='time', keep='last')

g = xr.open_dataset(grid)
lon = g.lon_rho
lat = g.lat_rho
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')

for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
    fig_temp, ax_temp = plt.subplots(figsize=(10, 5))
    fig_temp.suptitle(f'Temperature relative dispersion from the N+1 day spread in the {name} zone')  
    fig_salt, ax_salt = plt.subplots(figsize=(10, 5))
    fig_salt.suptitle(f'Salinity relative dispersion from the N+1 day spread in the {name} zone')   
    fig_mld, ax_mld = plt.subplots(figsize=(10, 5))
    fig_mld.suptitle(f'MLD relative dispersion from the N+1 day spread in the {name} zone')  
    fig_str, ax_str = plt.subplots(figsize=(10, 5))
    fig_str.suptitle(f'Wind stress relative dispersion from the N+1 day spread in the {name} zone')  
    
    print(f'├──Working on {name}.')
    box = (lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2)
    
    for (ens, std, mld_std, wstr_std) in zip(ensembles, [gls_std, str_std, ini_std], [gls_mld_std, str_mld_std, ini_mld_std], [gls_str_std, str_str_std, ini_str_std]):
        print(f'│   ├──Working on ensemble {ens}.')
        box_std = std.where(box.isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))).mean(dim=['xi_rho', 'eta_rho'], skipna=True)
        box_mld_std = mld_std.where(box.isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))).mean(dim=['xi_rho', 'eta_rho'], skipna=True)
        box_str_std = wstr_std.where(box.isel(xi_rho=slice(0, 431), eta_rho=slice(0, 501))).mean(dim=['xi_rho', 'eta_rho'], skipna=True)
        
        box_rel_std = box_std / box_std.isel(time=1)
        box_rel_mld_std = box_mld_std / box_mld_std.isel(time=1)
        box_rel_str_std = box_str_std / box_str_std.isel(time=1)
        
        ax_temp.semilogy(box_rel_std.time, box_rel_std.temp, color=ens_colors[ensembles.index(ens)], label=ens)
        ax_salt.semilogy(box_rel_std.time, box_rel_std.salt, color=ens_colors[ensembles.index(ens)], label=ens) 
        ax_mld.semilogy(box_rel_mld_std.time, box_rel_mld_std.mld, color=ens_colors[ensembles.index(ens)], label=ens)  
        ax_str.semilogy(box_rel_str_std.time, box_rel_str_std.windstr, color=ens_colors[ensembles.index(ens)], label=ens)   
    
        salt_dtime = box_rel_std.salt.isel(time=slice(1,100)).where(box_rel_std.salt>2, drop=True).isel(time=0).time
        temp_dtime = box_rel_std.temp.isel(time=slice(1,100)).where(box_rel_std.temp>2, drop=True).isel(time=0).time
        mld_dtime = box_rel_mld_std.mld.isel(time=slice(1,100)).where(box_rel_mld_std.mld>2, drop=True).isel(time=0).time
        str_dtime = box_rel_str_std.windstr.isel(time=slice(1,100)).where(box_rel_str_std.windstr>2, drop=True).isel(time=0).time
        
        salt_delta = pd.to_timedelta((salt_dtime - box_rel_std.time.isel(time=1)).data)
        temp_delta = pd.to_timedelta((temp_dtime - box_rel_std.time.isel(time=1)).data)
        mld_delta = pd.to_timedelta((mld_dtime - box_rel_mld_std.time.isel(time=1)).data)
        str_delta = pd.to_timedelta((str_dtime - box_rel_str_std.time.isel(time=1)).data)
        
        print(f'│   ├── Salinity doubling time: {salt_delta}.')
        print(f'│   ├── Temperature doubling time: {temp_delta}.')
        print(f'│   ├── MLD doubling time: {mld_delta}.')
        print(f'│   ├── Wind stress doubling time: {str_delta}.')
              
    ax_temp.set_xlim(np.datetime64('2016-12-10'), np.datetime64('2019-12-31'))
    ax_temp.set_ylim(0.5, 1000)
    ax_temp.set_ylabel('Temperature relative dispersion')
    ax_temp.tick_params("x", rotation=45)
    ax_temp.grid(linewidth=0.3)
    ax_temp.legend(loc='upper right')
    fig_temp.tight_layout()
              
    ax_salt.set_xlim(np.datetime64('2016-12-10'), np.datetime64('2019-12-31'))
    ax_salt.set_ylim(0.5, 1000)
    ax_salt.set_ylabel('Salinity relative dispersion')
    ax_salt.tick_params("x", rotation=45)
    ax_salt.grid(linewidth=0.3)
    ax_salt.legend(loc='upper right')
    fig_salt.tight_layout()
              
    ax_mld.set_xlim(np.datetime64('2016-12-10'), np.datetime64('2019-12-31'))
    ax_mld.set_ylim(0.5, 1000)
    ax_mld.set_ylabel('MLD relative dispersion')
    ax_mld.tick_params("x", rotation=45)
    ax_mld.grid(linewidth=0.3)
    ax_mld.legend(loc='upper right')
    fig_mld.tight_layout()
              
    ax_str.set_xlim(np.datetime64('2016-12-10'), np.datetime64('2019-12-31'))
    ax_str.set_ylim(0.05, 1000)
    ax_str.set_ylabel('MLD relative dispersion')
    ax_str.tick_params("x", rotation=45)
    ax_str.grid(linewidth=0.3)
    ax_str.legend(loc='upper right')
    fig_str.tight_layout()
    
    fig_temp.savefig(os.path.join(figures, f'all_ens_{name}_temp_relative_dispersion_self.png'), dpi=300, transparent=True)
    fig_salt.savefig(os.path.join(figures, f'all_ens_{name}_salt_relative_dispersion_self.png'), dpi=300, transparent=True)
    fig_mld.savefig(os.path.join(figures, f'all_ens_{name}_mld_relative_dispersion_self.png'), dpi=300, transparent=True)
    fig_str.savefig(os.path.join(figures, f'all_ens_{name}_str_relative_dispersion_self.png'), dpi=300, transparent=True)
    plt.show()
