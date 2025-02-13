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
from .utils import load_grid, load_data, save_figure, plot_map

def vel_vort_hel(data_path, start_date, end_date, figsize=(24, 8), cmap_velocity=cmcrameri.cm.oslo, cmap_vorticity=cmcrameri.cm.vik, cmap_helicity=cmcrameri.cm.vik):
    """
    Plot velocity, vorticity, and helicity data on a map.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Start time for the data slice.
    end_date : str
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
    u = u[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_3D_velocity(u, v, angle)
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
    ax.set_title(f"Velocity SWIO {start_date}", size=9)
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap_velocity.N)
    plot_map(ax, lon, lat, velocity, cmap_velocity, norm, levels, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

    # --- Vorticity Plot ---
    ax = axes[1]
    ax.set_title(f"Vorticity SWIO {start_date}", size=9)
    levels = np.linspace(-0.15, 0.15, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap_vorticity.N)
    plot_map(ax, lon, lat, vorticity * 3600, cmap_vorticity, norm, levels, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

    # --- Helicity Plot ---
    ax = axes[2]
    ax.set_title(f"Helicity SWIO {start_date}", size=9)
    levels = np.linspace(-0.5, 0.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap_helicity.N)
    plot_map(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap_helicity, norm, levels, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vel_vort_hel_{start_date}_{end_date}.png")
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
    u = u[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_3D_velocity(u, v, angle)
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(0, 2.5, 19)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, velocity, cmap, norm, levels, 'Velocity [$m.s^{-1}$]', msk, msk_inv, gridline_style)

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
    u = u[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')

    # Transform velocity components
    u_geo, v_geo = transform_3D_velocity(u, v, angle)

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
    plot_map(ax, lon, lat, vorticity * 3600, cmap, norm, levels, 'Vorticity [$h^{-1}$]', msk, msk_inv, gridline_style)

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
    u = u[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,-1,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    
    # Transform velocity components
    u_geo, v_geo = transform_3D_velocity(u, v, angle)
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
    plot_map(ax, lon, lat, helicity * 3600 ** 2 / 1000, cmap, norm, levels, 'Helicity [$km.h^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"helicity_{start_date}_{end_date}.png")
    plt.close(fig)
    
def eke(data_path, date, figsize=(8, 8), cmap=cmcrameri.cm.lapaz):
    """
    Plot EKE data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for EKE, by default cmcrameri.cm.lapaz
    """
    
    # Load grid data
    lon, lat, pm, pn, msk, msk_inv, angle = load_grid(is_Velocity=True)

    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan)
    v = v.where((v != fill_value), np.nan)
    w = w.where((w != fill_value), np.nan)
    print("NaN values added")
    
    # Moyenne annuelle
    u_yr = np.nanmean(u, axis = 0)
    v_yr = np.nanmean(v, axis = 0)
    w_yr = np.nanmean(w, axis = 0)
    print("Yearly mean calculated")
    
    u = u[:,:,:,:].sel(time=date)
    v = v[:,:,:,:].sel(time=date)
    w = w[:,:,:,:].sel(time=date)
    print("Date selected")
    
    # Vitesse turbulente
    ut = (u_yr - u).data
    vt = (v_yr - v).data
    wt = (w_yr - w).data
    print("Turbulent velocity calculated")
    
    angle_expand = angle[:,:].data.reshape(1, 1, angle.shape[0], angle.shape[1])
    
    ut_geo = ut[:,:,:-1,:] * np.cos(angle_expand[:,:,:-1,:-1]) - vt[:,:,:,:-1] * np.sin(angle_expand[:,:,:-1,:-1])
    vt_geo = ut[:,:,:-1,:] * np.sin(angle_expand[:,:,:-1,:-1]) + vt[:,:,:,:-1] * np.cos(angle_expand[:,:,:-1,:-1])
    wt_geo = wt[:,:,:-1,:-1]
    print("Turbulent velocity transformed")
    
    EKE = 1 / 2 * (ut_geo[0,:,:,:] ** 2 + vt_geo[0,:,:,:] ** 2 + wt_geo[0,:,:,:] ** 2)
    print("EKE calculated")
    
    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"EKE SWIO {date}", size=9)
    a = 1e-2
    b = 1
    c = 10
    levels = np.logspace(np.log10(a), np.log10(b), c * 2 - 1)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, EKE[-1,:,:], cmap, norm, levels, 'EKE [$m^2.s^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"eke_{date}.png")
    plt.close(fig)

def mke(data_path, start_date, end_date, figsize=(8, 8), cmap=cmcrameri.cm.lapaz):
    """
    Plot MKE data on a map for a specific date range.

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
    
    u = u[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    w = w[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    print("Data selected")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    w = w.where((w != fill_value), np.nan).data
    print("NaN values added")
    
    # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
    angle_expand = angle[:,:].data.reshape(1, angle.shape[0], angle.shape[1])
    
    u_geo = u[:,:-1,:] * np.cos(angle_expand[:,:-1,:-1]) - v[:,:,:-1] * np.sin(angle_expand[:,:-1,:-1])
    v_geo = u[:,:-1,:] * np.sin(angle_expand[:,:-1,:-1]) + v[:,:,:-1] * np.cos(angle_expand[:,:-1,:-1])
    w_geo = w[:,:-1,:-1]
    print("Velocity transformed")
    
    # Calculate EKE
    MKE = 1 / 2 * (u_geo ** 2 + v_geo ** 2 + w_geo ** 2)
    print("MKE calculated")

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"MKE SWIO {start_date} to {end_date}", size=9)
    a = 1e-2
    b = 1
    c = 10
    levels = np.logspace(np.log10(a), np.log10(b), c * 2 - 1)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon, lat, MKE[-1,:,:], cmap, norm, levels, 'MKE [$m^2.s^{-2}$]', msk, msk_inv, gridline_style)

    plt.tight_layout()
    save_figure(fig, f"mke_{start_date}_{end_date}.png")
    plt.close(fig)