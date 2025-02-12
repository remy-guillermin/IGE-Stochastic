"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage des données CROCO.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import xarray as xr
from .utils import load_grid, load_data, transform_velocity, save_figure, plot_data

def vel_vort_hel(data_path, start_time, end_time, figsize=(24, 8), cmap_velocity=cmcrameri.cm.oslo, cmap_vorticity=cmcrameri.cm.vik, cmap_helicity=cmcrameri.cm.vik):
    """
    Plot velocity, vorticity, and helicity data on a map.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_time : str
        Start time for the data slice.
    end_time : str
        End time for the data slice.
    figsize : tuple, optional
        Size of the figure, by default (24, 8)
    cmap_velocity : colormap, optional
        Colormap for velocity, by default cmcrameri.cm.oslo
    cmap_vorticity : colormap, optional
        Colormap for vorticity, by default cmcrameri.cm.vik
    cmap_helicity : colormap, optional
        Colormap for helicity, by default cmcrameri.cm.vik
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u, v, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=1) * pm
    du_dlat = np.gradient(u_geo, axis=0) * pn

    # Calculate vorticity and helicity
    vorticity = dv_dlon - du_dlat
    helicity = velocity * vorticity

    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    # --- Velocity Plot ---
    ax = axes[0]
    ax.set_title(f"Velocity SWIO {start_time}", size=9)
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap_velocity.N)
    plot_data(ax, lon, lat, velocity, cmap_velocity, norm, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    # --- Vorticity Plot ---
    ax = axes[1]
    ax.set_title(f"Vorticity SWIO {start_time}", size=9)
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap_vorticity.N)
    plot_data(ax, lon, lat, vorticity * 3600, cmap_vorticity, norm, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    # --- Helicity Plot ---
    ax = axes[2]
    ax.set_title(f"Helicity SWIO {start_time}", size=9)
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap_helicity.N)
    plot_data(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap_helicity, norm, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vel_vort_hel_{start_time}_{end_time}.png")
    plt.close(fig)


def velocity(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.oslo):
    """
    Plot velocity data on a map for a specific date range.

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
        Colormap for velocity, by default cmcrameri.cm.oslo
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u_mean, v_mean, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, velocity, cmap, norm, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"velocity_{start_date}_{end_date}.png")
    plt.close(fig)


def vorticity(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.vik):
    """
    Plot vorticity data on a map for a specific date range.

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
        Colormap for vorticity, by default cmcrameri.cm.vik
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u_mean, v_mean, angle)

    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=1) * pm
    du_dlat = np.gradient(u_geo, axis=0) * pn

    # Calculate vorticity
    vorticity = dv_dlon - du_dlat

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Vorticity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, vorticity * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vorticity_{start_date}_{end_date}.png")
    plt.close(fig)


def helicity(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.vik):
    """
    Plot helicity data on a map for a specific date range.

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
        Colormap for helicity, by default cmcrameri.cm.vik
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    
    # Transform velocity components
    u_geo, v_geo = transform_velocity(u_mean, v_mean, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=1) * pm
    du_dlat = np.gradient(u_geo, axis=0) * pn

    # Calculate vorticity and helicity
    vorticity = dv_dlon - du_dlat
    helicity = velocity * vorticity

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Helicity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap, norm, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"helicity_{start_date}_{end_date}.png")
    plt.close(fig)

def eke(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.lapaz):
    """
    Plot EKE data on a map for a specific date range.

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
        Colormap for EKE, by default cmcrameri.cm.lapaz
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
    u = u[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')
    w = w[:,-1,:,:].sel(time=slice(start_time, end_time)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u, v, angle)
    w_geo = w.data

    # Calculate EKE
    EKE = 1 / 2 * (u_geo ** 2 + v_geo ** 2 + w_geo ** 2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"EKE SWIO {start_date} to {end_date}", size=9)
    a = 1e-2
    b = 1
    c = 10
    levels = np.logspace(np.log10(a), np.log10(b), c * 2 - 1)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon[:-1,:-1], lat[:-1,:-1], EKE, cmap, norm, 'EKE [$m^2.s^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"eke_{start_date}_{end_date}.png")
    plt.close(fig)