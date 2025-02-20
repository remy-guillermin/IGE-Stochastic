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
from .utils import load_grid, load_data, save_figure, plot_map

def all_surf(data_path, start_date, end_date, figsize=(24, 6), sss_cmap=cmocean.cm.haline, ssh_cmap=cmcrameri.cm.roma_r, sst_cmap=cmocean.cm.thermal):
    """
    Plot SSS, SSH, and SST data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Start date for the data slice in 'YYYY-MM-DD' format.
    end_date : str
        End date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (24, 6)
    sss_cmap : colormap, optional
        Colormap for the SSS plot, by default cmocean.cm.haline
    ssh_cmap : colormap, optional
        Colormap for the SSH plot, by default cmcrameri.cm.roma_r
    sst_cmap : colormap, optional
        Colormap for the SST plot, by default cmocean.cm.thermal
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    salt, zeta, temp = load_data(data_path, ('salt', 'zeta', 'temp'))
    
    fill_value = 9.96921e+36
    salt = salt.where((salt != fill_value), np.nan)
    zeta = zeta.where((zeta != fill_value), np.nan)
    temp = temp.where((temp != fill_value), np.nan)
    
    salt = salt[:, -1, :, :].sel(time=slice(start_date, end_date)).mean(dim='time')
    zeta = zeta.sel(time=slice(start_date, end_date)).mean(dim='time')
    temp = temp[:, -1, :, :].sel(time=slice(start_date, end_date)).mean(dim='time')
    
    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    ax = axs[0]
    str_sss = r'$\overline{\overline{SSS}}$'
    ax.set_title(f"{str_sss} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(34, 36, 21)
    norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
    plot_map(ax, lon, lat, salt, sss_cmap, norm, levels, 'SSS [psu]', msk, msk_inv, gridline_style)
    
    ax = axs[1]
    str_ssh = r'$\overline{\overline{SSH}}$'
    ax.set_title(f"{str_ssh} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(0, 1, 21)
    norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
    plot_map(ax, lon, lat, zeta, ssh_cmap, norm, levels, 'SSH [m]', msk, msk_inv, gridline_style)
    
    ax = axs[2]
    str_sst = r'$\overline{\overline{SST}}$'
    ax.set_title(f"{str_sst} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(20, 30, 21)
    norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
    plot_map(ax, lon, lat, temp, sst_cmap, norm, levels, 'SST [°C]', msk, msk_inv, gridline_style)
    
    plt.tight_layout()
    save_figure(fig, f"sea_surface_{start_date}_{end_date}.png")
    plt.close(fig)
    
def all_annual_anomaly(data_path, figsize=(24, 6), sss_cmap=cmocean.cm.haline, ssh_cmap=cmcrameri.cm.roma_r, sst_cmap=cmocean.cm.thermal):
    """
    Plot annual SLA, STA, and SSA data on a map.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    figsize : tuple, optional
        Size of the figure, by default (24, 6)
    sss_cmap : colormap, optional
        Colormap for the SSS plot, by default cmocean.cm.haline
    ssh_cmap : colormap, optional
        Colormap for the SSH plot, by default cmcrameri.cm.roma_r
    sst_cmap : colormap, optional
        Colormap for the SST plot, by default cmocean.cm.thermal
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    salt, zeta, temp = load_data(data_path, ('salt', 'zeta', 'temp'))
    
    fill_value = 9.96921e+36
    salt = salt.where((salt != fill_value), np.nan)
    zeta = zeta.where((zeta != fill_value), np.nan)
    temp = temp.where((temp != fill_value), np.nan)
    
    salt_yr = salt[:, -1, :, :].mean(dim='time')
    zeta_yr = zeta.mean(dim='time')
    temp_yr = temp[:, -1, :, :].mean(dim='time')
    
    SLA = np.sqrt((zeta - zeta_yr)**2).mean(dim='time')
    STA = np.sqrt((temp - temp_yr)**2).mean(dim='time')
    SSA = np.sqrt((salt - salt_yr)**2).mean(dim='time')
    
    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    ax = axs[0]
    ax.set_title(r"$\overline{\overline{SSA}}$ SWIO", size=9)
    levels = np.linspace(0, 1, 11)
    norm = mpl.colors.BoundaryNorm(levels, sss_cmap.N)
    plot_map(ax, lon, lat, SSA[-1,:,:], sss_cmap, norm, levels, 'SSA [psu]', msk, msk_inv, gridline_style)
    
    ax = axs[1]
    ax.set_title(r"$\overline{\overline{SLA}}$ SWIO", size=9)
    levels = np.linspace(0, 0.2, 11)
    norm = mpl.colors.BoundaryNorm(levels, ssh_cmap.N)
    plot_map(ax, lon, lat, SLA, ssh_cmap, norm, levels, 'SLA [m]', msk, msk_inv, gridline_style)
    
    ax = axs[2]
    ax.set_title(r"$\overline{\overline{STA}}$ SWIO", size=9)
    levels = np.linspace(0, 5, 11)
    norm = mpl.colors.BoundaryNorm(levels, sst_cmap.N)
    plot_map(ax, lon, lat, STA[-1,:,:], sst_cmap, norm, levels, 'STA [°C]', msk, msk_inv, gridline_style)
    
    plt.tight_layout()
    save_figure(fig, f"sea_anomaly_annual_{data_path}.png")
    plt.close(fig)
    
def sss(data_path, start_date, end_date, figsize=(8, 8), cmap=cmocean.cm.haline):
    """
    Plot SSS data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Start date for the data slice in 'YYYY-MM-DD' format.
    end_date : str
        End date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for the plot, by default cmocean.cm.haline
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()

    # Load simulation data
    salt = load_data(data_path, ('salt',))[:, -1, :, :].sel(time=slice(start_date, end_date)).mean(dim='time')

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    str_sss = r'$\overline{\overline{SSS}}$'
    ax.set_title(f"{str_sss} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(34, 36, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, salt, cmap, norm, levels, 'SSS [psu]', msk, msk_inv, gridline_style)
    
    plt.tight_layout()
    save_figure(fig, f"sss_{start_date}_{end_date}.png")
    plt.close(fig)



def ssh(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.roma_r):
    """
    Plot SSH data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Start date for the data slice in 'YYYY-MM-DD' format.
    end_date : str
        End date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for the plot, by default cmcrameri.cm.roma_r
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()

    # Load simulation data
    zeta = load_data(data_path, ('zeta',)).sel(time=slice(start_date, end_date)).mean(dim='time')

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    str_ssh = r'$\overline{\overline{SSH}}$'
    ax.set_title(f"{str_ssh} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(0, 1, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, zeta, cmap, norm, levels, 'SSH [m]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"ssh_{start_date}_{end_date}.png")
    plt.close(fig)


def sst(data_path, start_date, end_date, figsize=(8, 8), cmap=cmocean.cm.thermal):
    """
    Plot SST data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Start date for the data slice in 'YYYY-MM-DD' format.
    end_date : str
        End date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for the plot, by default cmocean.cm.thermal
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()

    # Load simulation data
    temp = load_data(data_path, ('temp',))[:, -1, :, :].sel(time=slice(start_date, end_date)).mean(dim='time')

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}
    
    str_sst = r'$\overline{\overline{SST}}$'
    ax.set_title(f"{str_sst} SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(20, 30, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, temp, cmap, norm, levels, 'SST [°C]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"sst_{start_date}_{end_date}.png")
    plt.close(fig)