"""
Module for surface data visualization in CROCO.

This module provides functions to plot surface variables (e.g., SSS, SSH, SST) 
from CROCO simulation data or GLORYS datasets.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import os
from typing import Optional
from .utils import load_grid, load_data, load_GLORYS, save_figure, plot_map

def multiple_annual(data_path: str, 
                    variables: Optional[str] = 'all', 
                    start_date: Optional[str] = None, 
                    end_date: Optional[str] = None, 
                    figsize: Optional[tuple] = (24, 6), 
                    sss_cmap: Optional[mcolors.LinearSegmentedColormap] = cmocean.cm.haline, 
                    ssh_cmap: Optional[mcolors.ListedColormap] = cmcrameri.cm.roma_r, 
                    sst_cmap: Optional[mcolors.LinearSegmentedColormap] = cmocean.cm.thermal, 
                    grid_path: Optional[str] = None,
                    interactive: Optional[bool] = False):
    """
    Plot multiple surface variables (SSS, SSH, SST) on a map for a given date range 
    or compute their annual standard deviation.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    variables : list of str, optional
        Variables to plot. Options: 'sss', 'ssh', 'sst', 'annual_sd'. Default is 'all'.
    start_date : str, optional
        Start date for the data slice in 'YYYY-MM-DD' format. Default is None.
    end_date : str, optional
        End date for the data slice in 'YYYY-MM-DD' format. Default is None.
    figsize : tuple, optional
        Size of the figure in inches. Default is (24, 6).
    sss_cmap : colormap, optional
        Colormap for the SSS plot. Default is cmocean.cm.haline.
    ssh_cmap : colormap, optional
        Colormap for the SSH plot. Default is cmcrameri.cm.roma_r.
    sst_cmap : colormap, optional
        Colormap for the SST plot. Default is cmocean.cm.thermal.
    grid_path : str, optional
        Path to the grid data file. If None, the default grid is used. Default is None.
    interactive : bool, optional
        Whether to use an interactive backend for plotting. Default is False.
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
        start_date = str(np.datetime_as_string(time.data[0], 'D'))
        end_date = str(np.datetime_as_string(time.data[-1], 'D'))
    
    if 'annual_sd' in variables:
        salt_SD = salt.std(dim='time')
        temp_SD = temp.std(dim='time')
        zeta_SD = zeta.std(dim='time')
        
    salt = salt.mean(dim='time')
    zeta = zeta.mean(dim='time')
    temp = temp.mean(dim='time')
    
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    if 'sss' in variables and 'ssh' in variables and 'sst' in variables:
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

        ax = axs[0]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        levels = np.linspace(34, 36, 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt, sss_cmap, norm, r'$\overline{\overline{SSS}}$ [psu]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[1]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        levels = np.linspace(0, 1, 21)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta, ssh_cmap, norm, r'$\overline{\overline{SSH}}$ [m]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[2]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        levels = np.linspace(20, 30, 21)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp, sst_cmap, norm, r'$\overline{\overline{SST}}$ [°C]', msk, msk_inv, gridline_style, interactive=interactive)
        
        plt.tight_layout()
        save_figure(fig, f"surface_all_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
        
    if 'annual_sd':
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        ax = axs[0]
        ax.set_title("Annual Standard Deviation", size=9)
        levels = np.linspace(0, np.round(np.nanmax(salt_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt_SD, sss_cmap, norm, r'$SD_{SSS}$ [psu]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[1]
        levels = np.linspace(0, np.round(np.nanmax(zeta_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta_SD, ssh_cmap, norm, r'$SD_{SSH}$ [m]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[2]
        levels = np.linspace(0, np.round(np.nanmax(temp_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp_SD, sst_cmap, norm, r'$SD_{SST}$ [°C]', msk, msk_inv, gridline_style, interactive=interactive)

        plt.tight_layout()
        save_figure(fig, f"surface_SD_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
        
def GLORYS_annual(data_path: str,
                  variables: Optional[str] = 'all',
                  start_date: Optional[str] = None, 
                  end_date: Optional[str] = None, 
                  figsize: Optional[tuple] = (24, 6), 
                  sss_cmap: Optional[mcolors.LinearSegmentedColormap] = cmocean.cm.haline, 
                  ssh_cmap: Optional[mcolors.ListedColormap] = cmcrameri.cm.roma_r, 
                  sst_cmap: Optional[mcolors.LinearSegmentedColormap] = cmocean.cm.thermal,
                  interactive: Optional[bool] = False):
    """
    Plot GLORYS surface variables (SSS, SSH, SST) on a map for a given date range 
    or compute their annual standard deviation.

    Parameters
    ----------
    data_path : str
        Path to the GLORYS data file.
    variables : list of str, optional
        Variables to plot. Options: 'sss', 'ssh', 'sst', 'annual_sd'. Default is 'all'.
    start_date : str, optional
        Start date for the data slice in 'YYYY-MM-DD' format. Default is None.
    end_date : str, optional
        End date for the data slice in 'YYYY-MM-DD' format. Default is None.
    figsize : tuple, optional
        Size of the figure in inches. Default is (24, 6).
    sss_cmap : colormap, optional
        Colormap for the SSS plot. Default is cmocean.cm.haline.
    ssh_cmap : colormap, optional
        Colormap for the SSH plot. Default is cmcrameri.cm.roma_r.
    sst_cmap : colormap, optional
        Colormap for the SST plot. Default is cmocean.cm.thermal.
    interactive : bool, optional
        Whether to use an interactive backend for plotting. Default is False.
    """
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['sss', 'ssh', 'sst', 'annual_sd']
    
    salt, temp, zeta, _, _, lon, lat, msk, msk_inv = load_GLORYS(data_path)
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
        ax.set_title("Annual Standard Deviation")
        levels = np.linspace(0, np.round(np.nanmax(salt_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt_SD, sss_cmap, norm, r'$SD_{SSS}$ [psu]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[1]
        levels = np.linspace(0, np.round(np.nanmax(zeta_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta_SD, ssh_cmap, norm, r'$SD_{SSH}$ [m]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[2]
        levels = np.linspace(0, np.round(np.nanmax(temp_SD), 1), 21)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp_SD, sst_cmap, norm, r'$SD_{SST}$ [m]', msk, msk_inv, gridline_style, iteractive=interactive)

        plt.tight_layout()
        save_figure(fig, f"surface_SD_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
    
    if 'sss' in variables and 'ssh' in variables and 'sst' in variables:
        fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

        ax = axs[0]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        a = 33
        b = 36
        c = b - a
        levels = np.linspace(a, b, 4 * c + 1)
        norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
        plot_map(ax, lon, lat, salt, sss_cmap, norm, r'$\overline{\overline{SSS}}$ [psu]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[1]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        a = np.round(np.nanmin(zeta), 1)
        b = np.round(np.nanmax(zeta), 1)
        c = int((b - a) * 10)
        levels = np.linspace(a, b, c + 1)
        norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
        plot_map(ax, lon, lat, zeta, ssh_cmap, norm, r'$\overline{\overline{SSH}}$ [°C]', msk, msk_inv, gridline_style, interactive=interactive)
        
        ax = axs[2]
        ax.set_title(f"SWIO {start_date} to {end_date}")
        a = int(np.nanmin(temp))
        b = int(np.nanmax(temp))
        c = b - a
        levels = np.linspace(a, b, c + 1)
        norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
        plot_map(ax, lon, lat, temp, sst_cmap, norm, r'$\overline{\overline{SST}}$ [°C]', msk, msk_inv, gridline_style, interactive=interactive)
        
        plt.tight_layout()
        save_figure(fig, f"surface_all_{os.path.splitext(data_path)[0]}.png")
        plt.close(fig)
