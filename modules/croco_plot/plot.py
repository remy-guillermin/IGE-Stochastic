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
from .utils import load_grid, load_data, transform_velocity
import os
import subprocess  # Add this import


def plot_data(ax, lon, lat, data, cmap, norm, label, msk, msk_inv, gridline_style):
    """
    Helper function to plot data on a given axis.

    Parameters
    ----------
    ax : GeoAxes
        The axis to plot on.
    lon : ndarray
        Longitudes.
    lat : ndarray
        Latitudes.
    data : ndarray
        Data to plot.
    cmap : Colormap
        Colormap to use.
    norm : Normalize
        Normalization for the colormap.
    label : str
        Label for the colorbar.
    msk : ndarray
        Mask for contour.
    msk_inv : ndarray
        Inverse mask for contourf.
    gridline_style : dict
        Style for gridlines.
    """
    pcm = ax.pcolormesh(lon[:-1, :-1], lat[:-1, :-1], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
    ax.contourf(lon, lat, msk_inv, colors='lightgray')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}

    cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
    ticks = np.linspace(norm.vmin, norm.vmax, len(cb.get_ticks()))
    cb.set_ticks(ticks)
    cb.ax.set_yticklabels(np.round(ticks, 2), fontsize=8)


def save_figure(fig, filename):
    """
    Save the figure to the specified filename.

    Parameters
    ----------
    fig : Figure
        The figure to save.
    filename : str
        The path to save the figure.
    """
    output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, filename))
    print(f"Figure saved as {os.path.join(output_dir, filename)}.")


def open_figure(filename):
    """
    Open the saved figure using a terminal command without blocking the IPython session.

    Parameters
    ----------
    filename : str
        The name of the file to open.
    """
    output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
    file_path = os.path.join(output_dir, filename)
    if os.path.exists(file_path):
        subprocess.Popen(['eog', file_path])  # Use eog for Linux
    else:
        print(f"File {file_path} does not exist.")


def vel_vort_hel(data_path, start_time, end_time, figsize=(24, 8)):
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
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=slice(start_time, end_time))
    v = v.sel(time=slice(start_time, end_time))

    # Average over the selected time period
    u_mean = u.mean(dim='time')
    v_mean = v.mean(dim='time')

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
    fig, axes = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    # --- Velocity Plot ---
    ax = axes[0]
    ax.set_title(f"Velocity SWIO {start_time}", size=9)
    cmap = cmcrameri.cm.oslo
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, velocity, cmap, norm, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    # --- Vorticity Plot ---
    ax = axes[1]
    ax.set_title(f"Vorticity SWIO {start_time}", size=9)
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, vorticity * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    # --- Helicity Plot ---
    ax = axes[2]
    ax.set_title(f"Helicity SWIO {start_time}", size=9)
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap, norm, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vel_vort_hel_{start_time}_{end_time}.png")
    plt.close(fig)


def velocity(data_path, start_date, end_date, figsize=(8, 8)):
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
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=slice(start_date, end_date))
    v = v.sel(time=slice(start_date, end_date))

    # Average over the selected time period
    u_mean = u.mean(dim='time')
    v_mean = v.mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u_mean, v_mean, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {start_date} to {end_date}", size=9)
    cmap = cmcrameri.cm.oslo
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, velocity, cmap, norm, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"velocity_{start_date}_{end_date}.png")
    plt.close(fig)


def vorticity(data_path, start_date, end_date, figsize=(8, 8)):
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
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=slice(start_date, end_date))
    v = v.sel(time=slice(start_date, end_date))

    # Average over the selected time period
    u_mean = u.mean(dim='time')
    v_mean = v.mean(dim='time')

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
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, vorticity * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vorticity_{start_date}_{end_date}.png")
    plt.close(fig)


def helicity(data_path, start_date, end_date, figsize=(8, 8)):
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
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=slice(start_date, end_date))
    v = v.sel(time=slice(start_date, end_date))

    # Average over the selected time period
    u_mean = u.mean(dim='time')
    v_mean = v.mean(dim='time')
    
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
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap, norm, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"helicity_{start_date}_{end_date}.png")
    plt.close(fig)