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


def vel_vort_hel(data_path, start_time, end_time):
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
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})

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
    plt.show()


def velocity(data_path, date):
    """
    Plot velocity data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=date)
    v = v.sel(time=date)

    # Average over the selected time period
    u_mean = u.mean(dim='time')
    v_mean = v.mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_velocity(u_mean, v_mean, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {date}", size=9)
    cmap = cmcrameri.cm.oslo
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, velocity, cmap, norm, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    plt.show()


def vorticity(data_path, date):
    """
    Plot vorticity data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=date)
    v = v.sel(time=date)

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
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Vorticity SWIO {date}", size=9)
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, vorticity * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    plt.show()


def helicity(data_path, date):
    """
    Plot helicity data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    """
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    u = u.sel(time=date)
    v = v.sel(time=date)

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
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Helicity SWIO {date}", size=9)
    cmap = cmcrameri.cm.vik
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_data(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap, norm, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    plt.show()