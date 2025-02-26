"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage des données de surface CROCO.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import xarray as xr
import os
from .utils import load_grid, load_data, load_GLORYS, save_figure, plot_map

def multiple_annual(data_path, 
                    variables='all', 
                    start_date=None, 
                    end_date=None, 
                    figsize=(24, 6), 
                    sss_cmap=cmocean.cm.haline, 
                    ssh_cmap=cmcrameri.cm.roma_r, 
                    sst_cmap=cmocean.cm.thermal, 
                    grid_path=None):
    """
    Plot multiple surface variables (SSS, SSH, SST) on a map for a specific date range or annual standard deviation.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    variables : list of str, optional
        List of variables to plot, by default 'all'
        options : 'sss', 'ssh', 'sst', 'annual_sd'
    start_date : str, optional
        Start date for the data slice in 'YYYY-MM-DD' format.
    end_date : str, optional
        End date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (24, 6)
    sss_cmap : colormap, optional
        Colormap for the SSS plot, by default cmocean.cm.haline
    ssh_cmap : colormap, optional
        Colormap for the SSH plot, by default cmcrameri.cm.roma_r
    sst_cmap : colormap, optional
        Colormap for the SST plot, by default cmocean.cm.thermal
    
    cplot.time_series.multiple_time_series(['swio_avg_2017.nc', 'swio_avg_2018.nc', 'swio_avg_2019.nc', 'swio_avg_2020.nc', 'swio_avg_2021.nc', 'swio_avg_2022.nc', 'swio_avg_2023.nc'], variables='all')
    """
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['sss', 'ssh', 'sst', 'annual_sd']
    
    # Load grid data
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, _, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    salt, zeta, temp, time = load_data(data_path, ('salt', 'zeta', 'temp', 'time'))
    
    fill_value = 9.96921e+36
    salt = salt[:, -1, :, :]
    salt = salt.where((salt != fill_value), np.nan)
    zeta = zeta.where((zeta != fill_value), np.nan)
    temp = temp[:, -1, :, :]
    temp = temp.where((temp != fill_value), np.nan)
    
    if start_date and end_date:
        salt = salt.sel(time=slice(start_date, end_date))
        zeta = zeta.sel(time=slice(start_date, end_date))
        temp = temp.sel(time=slice(start_date, end_date))
    else:
        start_date = time.data[0]
        end_date = time.data[-1]
    
    if 'annual_sd' in variables:
        salt_yr = salt.mean(dim='time')
        zeta_yr = zeta.mean(dim='time')
        temp_yr = temp.mean(dim='time')
        
        zeta_SD = np.sqrt((zeta - zeta_yr)**2).mean(dim='time')
        temp_SD = np.sqrt((temp - temp_yr)**2).mean(dim='time')
        salt_SD = np.sqrt((salt - salt_yr)**2).mean(dim='time')
        
    salt = salt.mean(dim='time')
    zeta = zeta.mean(dim='time')
    temp = temp.mean(dim='time')
    
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    if 'sss' in variables and 'ssh' in variables and 'sst' in variables:
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

        ax = axs[0]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        levels = np.linspace(34, 36, 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt, sss_cmap, norm, r'$\overline{\overline{SSS}}$ [psu]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[1]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        levels = np.linspace(0, 1, 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, zeta, ssh_cmap, norm, r'$\overline{\overline{SSH}}$ [m]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[2]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        levels = np.linspace(20, 30, 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, temp, sst_cmap, norm, r'$\overline{\overline{SST}}$ [°C]', msk, msk_inv, gridline_style, levels = levels)
        
        plt.tight_layout()
        save_figure(fig, f"surface_all_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
        
    if 'annual_sd':
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        ax = axs[0]
        ax.set_title("Annual Standard Deviation", size=9)
        levels = np.linspace(0, np.round(np.nanmax(salt_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt_SD, sss_cmap, norm, r'$SD_{SSS}$ [psu]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[1]
        levels = np.linspace(0, np.round(np.nanmax(zeta_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta_SD, ssh_cmap, norm, r'$SD_{SSH}$ [m]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[2]
        levels = np.linspace(0, np.round(np.nanmax(temp_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp_SD, sst_cmap, norm, r'$SD_{SST}$ [°C]', msk, msk_inv, gridline_style, levels = levels)

        plt.tight_layout()
        save_figure(fig, f"surface_SD_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
        
def GLORYS_annual(
    data_path: str,
    variables: list = 'all',
    start_date: str = None, 
    end_date: str = None, 
    figsize=(24, 6), 
    sss_cmap=cmocean.cm.haline,
    ssh_cmap=cmcrameri.cm.roma_r, 
    sst_cmap=cmocean.cm.thermal, 
):
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['sss', 'ssh', 'sst', 'annual_sd']
    
    salt, temp, zeta, u, v, lon, lat, msk, msk_inv = load_GLORYS(data_path)
    salt = salt[:,0,:,:]
    zeta = zeta
    temp = temp[:,0,:,:]
    
    if 'annual_sd' in variables:
        salt_yr = salt.mean(dim='time')
        zeta_yr = zeta.mean(dim='time')
        temp_yr = temp.mean(dim='time')
        
        temp_SD = np.sqrt((temp - temp_yr)**2).mean(dim='time')
        zeta_SD = np.sqrt((zeta - zeta_yr)**2).mean(dim='time')
        salt_SD = np.sqrt((salt - salt_yr)**2).mean(dim='time')
        
    salt = salt.mean(dim='time')
    zeta = zeta.mean(dim='time')
    temp = temp.mean(dim='time')
    
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    if 'annual_sd' in variables:
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        ax = axs[0]
        ax.set_title("Annual Standard Deviation", size=9)
        levels = np.linspace(0, np.round(np.nanmax(salt_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt_SD, sss_cmap, norm, r'$SD_{SSS}$ [psu]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[1]
        levels = np.linspace(0, np.round(np.nanmax(zeta_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta_SD, ssh_cmap, norm, r'$SD_{SSH}$ [m]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[2]
        levels = np.linspace(0, np.round(np.nanmax(temp_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp_SD, sst_cmap, norm, r'$SD_{SST}$ [m]', msk, msk_inv, gridline_style, levels = levels)

        plt.tight_layout()
        save_figure(fig, f"surface_SD_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
    
    if 'sss' in variables and 'ssh' in variables and 'sst' in variables:
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

        ax = axs[0]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        a = 33
        b = 36
        c = b - a
        levels = np.linspace(a, b, 4 * c + 1)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt, sss_cmap, norm, r'$\overline{\overline{SSS}}$ [psu]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[1]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        a = np.round(np.nanmin(zeta), 1)
        b = np.round(np.nanmax(zeta), 1)
        c = int((b - a) * 10)
        levels = np.linspace(a, b, c + 1)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta, ssh_cmap, norm, r'$\overline{\overline{SSH}}$ [°C]', msk, msk_inv, gridline_style, levels = levels)
        
        ax = axs[2]
        ax.set_title(f"SWIO {start_date} to {end_date}", size=9)
        a = int(np.nanmin(temp))
        b = int(np.nanmax(temp))
        c = b - a
        levels = np.linspace(a, b, c + 1)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp, sst_cmap, norm, r'$\overline{\overline{SST}}$ [°C]', msk, msk_inv, gridline_style, levels = levels)
        
        plt.tight_layout()
        save_figure(fig, f"surface_all_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
