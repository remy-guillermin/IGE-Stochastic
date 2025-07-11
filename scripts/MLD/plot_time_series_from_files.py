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

lat_index = [430, 280]

ensembles = ['GLS', 'STR', 'INI']
ens_colors = ['crimson', 'royalblue', 'black']

box_colors = ['sandybrown', 'slateblue', 'lightseagreen', 'red']
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']



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



rms_list = {ensemble : {} for ensemble in ensembles}
deviation_list = {name : {} for name in names}
mean_list = {name : {} for name in names}
member_mean_list = {name : {} for name in names}

for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
    print(f'├──Working on box {name}.')
    # Box limits
    region_mask = (lon > lon1) & (lon < lon2) & (lat > lat1) & (lat < lat2)
    for ds_mean, ds_std, ens in zip([gls_mean, str_mean, ini_mean], [gls_std, str_std, ini_std], ensembles):
        print(f'│   ├──Working on {ens} ensemble.')
        ensemble_mean2d = ds_mean.mld.where(region_mask, drop=True)
        ensemble_std2d = ds_std.mld.where(region_mask, drop=True)
        rms_list[ens][name] = np.sqrt((ensemble_std2d ** 2).mean(dim=['eta_rho', 'xi_rho'], skipna=True))
        deviation_list[name][ens] = ensemble_std2d.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
        mean_list[name][ens] = ensemble_mean2d.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
        
    print('│   └──Done.')
print('Finished.')



for i, (box, box_deviation_list) in enumerate(deviation_list.items()):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f'MLD ensemble standard deviation in {box} box')
    for j, (ens, box_deviation) in enumerate(box_deviation_list.items()):
        ax.fill_between(
            box_deviation.time,
            - box_deviation,
            box_deviation,
            alpha=0.5,
            color=ens_colors[j],
            label=ens
        )
        
    ax.set_ylim(-40, 40)
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
    ax.set_ylabel('MLD deviation (m)')
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.5)
    
    fig.autofmt_xdate()
    fig.savefig(os.path.join(figures, f'mld_std_{box}.png'), dpi=300, transparent=True)
    plt.close()



for i, (box, box_mean_list) in enumerate(mean_list.items()):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f'MLD ensemble mean in {box} box')
    for j, (ens, box_mean) in enumerate(box_mean_list.items()):
        ax.plot(box_mean.time, box_mean, color=ens_colors[j], label=ens, linewidth=0.5)
        

    ax.set_ylabel('MLD mean (m)')
    ax.set_ylim(-120, 0)
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
    ax.legend(loc='upper right')
    
    fig.autofmt_xdate()
    fig.savefig(os.path.join(figures, f'mld_mean_{box}.png'), dpi=300, transparent=True)
    plt.close()
