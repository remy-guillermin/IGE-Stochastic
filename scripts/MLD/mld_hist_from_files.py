# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
from scipy.stats import norm
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



i = 200 #410
HIST_BINS = np.linspace(-150, 0, 20)
date = str(str_mean.isel(time=i).time.values.astype('datetime64[D]'))

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

for ax, (lon1, lon2, lat1, lat2), name in zip(axs.flatten(), boxes, names):
    mask = (lon > lon1) & (lon < lon2) & (lat > lat1) & (lat < lat2)
    
    ax.set_title(name)
    
    data_ini = np.array([])
    for ds in ini_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_ini = np.concatenate((data_ini, data_ds))
    
    data_str = np.array([])
    for ds in str_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_str = np.concatenate((data_str, data_ds))
    
    data_gls = np.array([])
    for ds in gls_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_gls = np.concatenate((data_gls, data_ds))
    
    n, bins, bar_container = ax.hist([data_ini, data_str, data_gls], HIST_BINS, density=True, color=colors, label=[f'INI', f'STR', f'GLS'], lw=1, alpha=0.5, edgecolor='black')
    
    xmin, xmax = ax.get_xlim()
    
    # x = np.linspace(xmin, xmax, 200)
    # p_ini = norm.pdf(x, np.nanmean(data_ini), np.nanstd(data_ini))
    # p_str = norm.pdf(x, np.nanmean(data_str), np.nanstd(data_str))
    # p_gls = norm.pdf(x, np.nanmean(data_gls), np.nanstd(data_gls))
    
    # ax.plot(x, p_ini, colors[0])
    # ax.plot(x, p_str, colors[1])
    # ax.plot(x, p_gls, colors[2])
    
    ax.set_xlabel('MLD (m)')
    ax.set_ylabel('Density')
    
    ax.set_ylim(0, 0.09)
    
    ini_text = AnchoredText(f'INI\nmean = {np.nanmean(data_ini):.1f}\nstd = {np.nanstd(data_ini):.1f}', 2, frameon=True)
    ini_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    ini_text.patch.set_alpha(0.5)
    ini_text.patch.set_edgecolor('black')
    ini_text.patch.set_linewidth(1.0)
    ax.add_artist(ini_text)
    
    str_text = AnchoredText(f'STR\nmean = {np.nanmean(data_str):.1f}\nstd = {np.nanstd(data_str):.1f}', 9, frameon=True)
    str_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    str_text.patch.set_facecolor(colors[1])
    str_text.patch.set_alpha(0.5)
    str_text.patch.set_edgecolor('black')
    str_text.patch.set_linewidth(1.0)
    ax.add_artist(str_text)
    
    gls_text = AnchoredText(f'GLS\nmean = {np.nanmean(data_gls):.1f}\nstd = {np.nanstd(data_gls):.1f}', 1, frameon=True)
    gls_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    gls_text.patch.set_facecolor(colors[2])
    gls_text.patch.set_alpha(0.5)
    gls_text.patch.set_edgecolor('black')
    gls_text.patch.set_linewidth(1.0)
    ax.add_artist(gls_text)
    
    print(f'Stats for {name}:')
    print(f'INI: mean = {np.nanmean(data_ini):.2f}, std = {np.nanstd(data_ini):.2f}')
    print(f'STR: mean = {np.nanmean(data_str):.2f}, std = {np.nanstd(data_str):.2f}')
    print(f'GLS: mean = {np.nanmean(data_gls):.2f}, std = {np.nanstd(data_gls):.2f}')
    
fig.suptitle(date)

fig.tight_layout()
fig.savefig(f'{figures}/hist_ensemble_{date}.png', dpi=300, transparent=True)
plt.show()



i = 410
HIST_BINS = np.linspace(-150, 0, 20)
date = str(str_mean.isel(time=i).time.values.astype('datetime64[D]'))

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

for ax, (lon1, lon2, lat1, lat2), name in zip(axs.flatten(), boxes, names):
    mask = (lon > lon1) & (lon < lon2) & (lat > lat1) & (lat < lat2)
    
    ax.set_title(name)
    
    data_ini = np.array([])
    for ds in ini_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_ini = np.concatenate((data_ini, data_ds))
    
    data_str = np.array([])
    for ds in str_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_str = np.concatenate((data_str, data_ds))
    
    data_gls = np.array([])
    for ds in gls_ds:
        data_ds = ds.sel(time=date).isel(time=0).mld
        data_ds = data_ds.where(mask, drop=True).values.flatten()
        data_ds = data_ds[~np.isnan(data_ds)]
        data_gls = np.concatenate((data_gls, data_ds))
    
    n, bins, bar_container = ax.hist([data_ini, data_str, data_gls], HIST_BINS, density=True, color=colors, lw=1, alpha=0.5, edgecolor='black')
    
    xmin, xmax = ax.get_xlim()
    
    # x = np.linspace(xmin, xmax, 200)
    # p_ini = norm.pdf(x, np.nanmean(data_ini), np.nanstd(data_ini))
    # p_str = norm.pdf(x, np.nanmean(data_str), np.nanstd(data_str))
    # p_gls = norm.pdf(x, np.nanmean(data_gls), np.nanstd(data_gls))
    
    # pdf_ini = ax.plot(x, p_ini, colors[0])
    # pdf_str = ax.plot(x, p_str, colors[1])
    # pdf_gls = ax.plot(x, p_gls, colors[2])
    
    ax.set_xlabel('MLD (m)')
    ax.set_ylabel('Density')
    
    ax.set_ylim(0, 0.09)
    
    ini_text = AnchoredText(f'INI\nmean = {np.nanmean(data_ini):.1f}\nstd = {np.nanstd(data_ini):.1f}', 2, frameon=True)
    ini_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    ini_text.patch.set_alpha(0.5)
    ini_text.patch.set_edgecolor('black')
    ini_text.patch.set_linewidth(1.0)
    ax.add_artist(ini_text)
    
    str_text = AnchoredText(f'STR\nmean = {np.nanmean(data_str):.1f}\nstd = {np.nanstd(data_str):.1f}', 9, frameon=True)
    str_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    str_text.patch.set_facecolor(colors[1])
    str_text.patch.set_alpha(0.5)
    str_text.patch.set_edgecolor('black')
    str_text.patch.set_linewidth(1.0)
    ax.add_artist(str_text)
    
    gls_text = AnchoredText(f'GLS\nmean = {np.nanmean(data_gls):.1f}\nstd = {np.nanstd(data_gls):.1f}', 1, frameon=True)
    gls_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    gls_text.patch.set_facecolor(colors[2])
    gls_text.patch.set_alpha(0.5)
    gls_text.patch.set_edgecolor('black')
    gls_text.patch.set_linewidth(1.0)
    ax.add_artist(gls_text)
    
    print(f'Stats for {name}:')
    print(f'INI: mean = {np.nanmean(data_ini):.2f}, std = {np.nanstd(data_ini):.2f}')
    print(f'STR: mean = {np.nanmean(data_str):.2f}, std = {np.nanstd(data_str):.2f}')
    print(f'GLS: mean = {np.nanmean(data_gls):.2f}, std = {np.nanstd(data_gls):.2f}')
    
fig.suptitle(date)

fig.tight_layout()
fig.savefig(f'{figures}/hist_ensemble_{date}.png', dpi=300, transparent=True)
plt.show()










