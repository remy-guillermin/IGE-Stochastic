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

ensembles = ['GLS', 'STR', 'INI']
ens_colors = ['crimson', 'royalblue', 'black']
ens_linestyles = ['solid', 'dashed', 'dotted']

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



ini_ds = [xr.open_dataset(f).sortby('time') for f in ini_files]
ini_mean = xr.open_dataset(ini_mean_file).sortby('time')
ini_std = xr.open_dataset(ini_std_file).sortby('time')

str_ds = [xr.open_dataset(f).sortby('time') for f in str_files]
str_mean = xr.open_dataset(str_mean_file).sortby('time')
str_std = xr.open_dataset(str_std_file).sortby('time')

gls_ds = [xr.open_dataset(f).sortby('time') for f in gls_files]
gls_mean = xr.open_dataset(gls_mean_file).sortby('time')
gls_std = xr.open_dataset(gls_std_file).sortby('time')

g = xr.open_dataset(grid)
lon = g.lon_rho
lat = g.lat_rho
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')




for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f'MLD deviation from the ensemble mean in the {name} zone')  
    
    print(f'├──Working on {name}.')
    box = (lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2)
    for (ens, mean, std) in zip(ensembles, [gls_mean,str_mean, ini_mean], [gls_std, str_std, ini_std]):
        print(f'│   ├──Working on ensemble {ens}.')
        box_mean = mean.where(box).mean(dim=['xi_rho', 'eta_rho'], skipna=True)
        box_std = std.where(box).mean(dim=['xi_rho', 'eta_rho'], skipna=True)
        ax.fill_between(box_std.time, -box_std.mld, box_std.mld, alpha=0.5, color=ens_colors[ensembles.index(ens)], label=ens)
              
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
    # ax.set_ylim(-0.75, 0.75)
    ax.set_ylabel('MLD deviation [m]')
    ax.tick_params("x", rotation=45)
    ax.grid(linewidth=0.3)
    ax.legend(loc='upper right')
    fig.tight_layout()
    
    fig.savefig(os.path.join(figures, f'all_ens_{name}_mld_deviation.png'), dpi=300, transparent=True)
    plt.show()
    






