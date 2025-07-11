# LIBRAIRIES
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import cmocean
import cartopy.crs as ccrs
import glob
import os
from functools import partial
import scipy.stats as sts


figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio.nc'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

lat_index = [430, 280]

ensembles = ['INI', 'STR', 'GLS']
ens_colors = ['black', 'royalblue', 'crimson']
ens_linestyles = ['solid', 'dashed', 'dotted']


ini_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'mean*.nc'))[0]
print(f'ini mld mean file found: {ini_mean_file}')
ini_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_ini', 'std*.nc'))[0]
print(f'ini mld std file found: {ini_std_file}')

str_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'mean*.nc'))[0]
print(f'str mld mean file found: {str_mean_file}')
str_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_str', 'std*.nc'))[0]
print(f'str mld std file found: {str_std_file}')

gls_mean_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'mean*.nc'))[0]
print(f'gls mld mean file found: {gls_mean_file}')
gls_std_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD/run_swio2_stoens30_gls', 'std*.nc'))[0]
print(f'gls mld std file found: {gls_std_file}')

ini_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_ini', 'mean*str.nc'))[0]
print(f'ini str mean file found: {ini_mean_str_file}')
ini_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_ini', 'std*str.nc'))[0]
print(f'ini str std file found: {ini_std_str_file}')

str_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_str', 'mean*str.nc'))[0]
print(f'str str mean file found: {str_mean_str_file}')
str_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_str', 'std*str.nc'))[0]
print(f'str str std file found: {str_std_str_file}')

gls_mean_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_gls', 'mean*str.nc'))[0]
print(f'gls str mean file found: {gls_mean_str_file}')
gls_std_str_file = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_gls', 'std*str.nc'))[0]
print(f'gls str std file found: {gls_std_str_file}')

ini_mean = xr.open_dataset(ini_mean_file).sortby('time')
ini_std = xr.open_dataset(ini_std_file).sortby('time')

ini_mean = ini_mean.drop_duplicates(dim='time', keep='last')
ini_std = ini_std.drop_duplicates(dim='time', keep='last')

str_mean = xr.open_dataset(str_mean_file).sortby('time')
str_std = xr.open_dataset(str_std_file).sortby('time')

str_mean = str_mean.drop_duplicates(dim='time', keep='last')
str_std = str_std.drop_duplicates(dim='time', keep='last')

gls_mean = xr.open_dataset(gls_mean_file).sortby('time')
gls_std = xr.open_dataset(gls_std_file).sortby('time')

gls_mean = gls_mean.drop_duplicates(dim='time', keep='last')
gls_std = gls_std.drop_duplicates(dim='time', keep='last')

ini_str_mean = xr.open_dataset(ini_mean_str_file).sortby('time')
ini_str_std = xr.open_dataset(ini_std_str_file).sortby('time')

ini_str_mean = ini_str_mean.drop_duplicates(dim='time', keep='last')
ini_str_std = ini_str_std.drop_duplicates(dim='time', keep='last')

str_str_mean = xr.open_dataset(str_mean_str_file).sortby('time')
str_str_std = xr.open_dataset(str_std_str_file).sortby('time')

str_str_mean = str_str_mean.drop_duplicates(dim='time', keep='last')
str_str_std = str_str_std.drop_duplicates(dim='time', keep='last')

gls_str_mean = xr.open_dataset(gls_mean_str_file).sortby('time')
gls_str_std = xr.open_dataset(gls_std_str_file).sortby('time')

gls_str_mean = gls_str_mean.drop_duplicates(dim='time', keep='last')
gls_str_std = gls_str_std.drop_duplicates(dim='time', keep='last')

g = xr.open_dataset(grid)
lon = g.lon_rho
lat = g.lat_rho
h = g.h
print(f'grid file opened: {grid}')
depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')

fig, axes = plt.subplots(2, 4, figsize=(21, 10))

for ax in axes.flatten():
    ax.grid(linewidth=0.3)

for ax in axes[:,0]:
    ax.set_ylabel('PDF')
    
HIST_BINS = np.linspace(-0, 0.03, 25)
for ax, index in zip(axes[0,:], [1, 10, 90, 900]):
    ax.set_xlabel(r'Wind stress (N m$^{-2}$)')
    ax.set_title(f'{index} days after start', fontsize=12)

    datas = [
        ini_str_std.isel(time=index).windstr.data.flatten(), 
        str_str_std.isel(time=index).windstr.data.flatten(), 
        gls_str_std.isel(time=index).windstr.data.flatten()
    ]
    
    ax.hist(
        datas, 
        HIST_BINS, 
        weights=[np.ones_like(data) / np.sum(~np.isnan(data)) for data in datas],
        color=ens_colors,  
        lw=1, alpha=0.5, 
        edgecolor='black'
    ) 
    
    ini_kde = sts.gaussian_kde(datas[0][~np.isnan(datas[0])])
    ax.plot(x, ini_kde.pdf(x), c=ens_colors[0], lw=2, label='INI')
    
    str_kde = sts.gaussian_kde(datas[1][~np.isnan(datas[1])])
    ax.plot(x, str_kde.pdf(x), c=ens_colors[1], lw=2, label='STR')
    
    gls_kde = sts.gaussian_kde(datas[2][~np.isnan(datas[2])])
    ax.plot(x, gls_kde.pdf(x), c=ens_colors[2], lw=2, label='GLS')   
    
    ini_text = AnchoredText(f'INI\nmean = {np.nanmean(datas[0])*1000:.1f}\nstd = {np.nanstd(datas[0])*1000:.1f}', 1, frameon=True)
    ini_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    ini_text.patch.set_alpha(0.5)
    ini_text.patch.set_edgecolor('black')
    ini_text.patch.set_linewidth(1.0)
    ax.add_artist(ini_text)
    
    str_text = AnchoredText(f'STR\nmean = {np.nanmean(datas[1])*1000:.1f}\nstd = {np.nanstd(datas[1])*1000:.1f}', 9, frameon=True)
    str_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    str_text.patch.set_facecolor(ens_colors[1])
    str_text.patch.set_alpha(0.5)
    str_text.patch.set_edgecolor('black')
    str_text.patch.set_linewidth(1.0)
    ax.add_artist(str_text)
    
    gls_text = AnchoredText(f'GLS\nmean = {np.nanmean(datas[2])*1000:.1f}\nstd = {np.nanstd(datas[2])*1000:.1f}', 2, frameon=True)
    gls_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    gls_text.patch.set_facecolor(ens_colors[2])
    gls_text.patch.set_alpha(0.5)
    gls_text.patch.set_edgecolor('black')
    gls_text.patch.set_linewidth
    ax.add_artist(gls_text)
        
HIST_BINS = np.linspace(0, 20, 25)
x = np.linspace(0, 20, 100)
for ax, index in zip(axes[1,:], [1, 10, 90, 900]):
    ax.set_xlabel('MLD (m)')

    datas = [
        ini_std.isel(time=index).mld.data.flatten(), 
        str_std.isel(time=index).mld.data.flatten(), 
        gls_std.isel(time=index).mld.data.flatten()
    ]
    
    ax.hist(
        datas, 
        HIST_BINS, 
        weights=[np.ones_like(data) / np.sum(~np.isnan(data)) for data in datas],
        color=ens_colors,
        lw=1, alpha=0.5, 
        edgecolor='black'
    )
    
    ini_kde = sts.gaussian_kde(datas[0][~np.isnan(datas[0])])
    ax.plot(x, ini_kde.pdf(x), c=ens_colors[0], lw=2, label='INI')
    
    str_kde = sts.gaussian_kde(datas[1][~np.isnan(datas[1])])
    ax.plot(x, str_kde.pdf(x), c=ens_colors[1], lw=2, label='STR')
    
    gls_kde = sts.gaussian_kde(datas[2][~np.isnan(datas[2])])
    ax.plot(x, gls_kde.pdf(x), c=ens_colors[2], lw=2, label='GLS')
    
    ini_text = AnchoredText(f'INI\nmean = {np.nanmean(datas[0]):.1f}\nstd = {np.nanstd(datas[0]):.1f}', 1, frameon=True)
    ini_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    ini_text.patch.set_alpha(0.5)
    ini_text.patch.set_edgecolor('black')
    ini_text.patch.set_linewidth(1.0)
    ax.add_artist(ini_text)
    
    str_text = AnchoredText(f'STR\nmean = {np.nanmean(datas[1]):.1f}\nstd = {np.nanstd(datas[1]):.1f}', 9, frameon=True)
    str_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    str_text.patch.set_facecolor(ens_colors[1])
    str_text.patch.set_alpha(0.5)
    str_text.patch.set_edgecolor('black')
    str_text.patch.set_linewidth(1.0)
    ax.add_artist(str_text)
    
    gls_text = AnchoredText(f'GLS\nmean = {np.nanmean(datas[2]):.1f}\nstd = {np.nanstd(datas[2]):.1f}', 2, frameon=True)
    gls_text.patch.set_boxstyle("round,pad=0.3,rounding_size=0.2")
    gls_text.patch.set_facecolor(ens_colors[2])
    gls_text.patch.set_alpha(0.5)
    gls_text.patch.set_edgecolor('black')
    gls_text.patch.set_linewidth
    ax.add_artist(gls_text)
    
for ax in axes[0,:]:
    ax.set_xlim(0, 0.025)
    ax.set_ylim(0.0, 1.5)
    
for ax in axes[1,:]:
    ax.set_ylim(0.0, 1.5)

fig.tight_layout()
fig.savefig(os.path.join(figures, 'windstr_mld_histograms.png'), dpi=300, transparent=True)
plt.show()  