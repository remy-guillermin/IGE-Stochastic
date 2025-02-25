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
from .utils import load_grid, load_data, save_figure, plot_map

def velocity(data_path, 
             start_date, 
             end_date, 
             figsize=(8, 8), 
             cmap=cmcrameri.cm.oslo, 
             grid_path=None):
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
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    u = u[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    print("Data sliced")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    print("NaN values added")
    
    # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
    angle_expand = angle[:,:].data.reshape(1, angle.shape[0], angle.shape[1])
    
    u_geo = u[:,:-1,:] * np.cos(angle_expand[:,:-1,:-1]) - v[:,:,:-1] * np.sin(angle_expand[:,:-1,:-1])
    v_geo = u[:,:-1,:] * np.sin(angle_expand[:,:-1,:-1]) + v[:,:,:-1] * np.cos(angle_expand[:,:-1,:-1])
    print("Velocity transformed")
    
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(0, 2.5, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity[-1,:,:], cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style, levels=levels)

    plt.tight_layout()
    save_figure(fig, f"velocity_{start_date}_{end_date}.png")
    plt.close(fig)


def vorticity(data_path, 
              start_date, 
              end_date, 
              figsize=(8, 8), 
              cmap=cmcrameri.cm.vik, 
              grid_path=None):
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
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    u = u[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    v = v[:,:,:,:].sel(time=slice(start_date, end_date)).mean(dim='time')
    print("Data selected")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    print("NaN values added")
    
    # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
    angle_expand = angle[:,:].data.reshape(1, angle.shape[0], angle.shape[1])
    
    u_geo = u[:,:-1,:] * np.cos(angle_expand[:,:-1,:-1]) - v[:,:,:-1] * np.sin(angle_expand[:,:-1,:-1])
    v_geo = u[:,:-1,:] * np.sin(angle_expand[:,:-1,:-1]) + v[:,:,:-1] * np.cos(angle_expand[:,:-1,:-1])
    print("Velocity transformed")

    pm_expand = pm.data.reshape(1, pm.shape[0], pm.shape[1])
    pn_expand = pn.data.reshape(1, pn.shape[0], pn.shape[1])
    
    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=2) * pm_expand[:,:-1,:-1]
    du_dlat = np.gradient(u_geo, axis=1) * pn_expand[:,:-1,:-1]

    # Calculate vorticity and helicity
    vorticity = dv_dlon - du_dlat

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Vorticity SWIO {start_date} to {end_date}", size=9)
    levels = np.linspace(-0.15, 0.15, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], vorticity[-1,:,:] * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style, levels=levels)

    plt.tight_layout()
    save_figure(fig, f"vorticity_{start_date}_{end_date}.png")
    plt.close(fig)
    
def eke(data_path, 
        date, 
        figsize=(8, 8), 
        cmap=cmcrameri.cm.lapaz, 
        grid_path=None):
    """
    Plot EKE data on a map for a specific date.

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
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, h = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, h = load_grid()

    # Load simulation data
    u, v, w, s_rho = load_data(data_path, ('u', 'v', 'w', 's_rho'))
    
    depth = h * s_rho # Profondeur
    cell_volume = -(depth[:,:,-1] * 1 / pn * 1 / pm).data # Volumes des cellules de surface
    surface_volume = - np.sum(depth[:,:,-1].data[h.data != 50] * 1 / pn.data[h.data != 50] * 1 / pm.data[h.data != 50])
    domain_volume = np.sum(h.data[h.data != 50] * 1/ pn.data[h.data != 50] * 1 / pm.data[h.data != 50]) # Volume du domaine
    depth, h ,s_rho, pm, pn = None, None, None, None, None 
    
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
    
    EKE = 1 / 2 * (ut_geo[0,-1,:,:] ** 2 + vt_geo[0,-1,:,:] ** 2 + wt_geo[0,-1,:,:] ** 2) * cell_volume[:-1,:-1] / surface_volume
    print("EKE calculated")
    
    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"EKE SWIO {date}", size=9)
    a = int(np.log10(np.nanmax(EKE)))
    b = a-2
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], EKE, cmap, norm, 'EKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"eke_{date}.png")
    plt.close(fig)

def mke(data_path, 
        start_date, 
        end_date, 
        figsize=(8, 8), 
        cmap=cmcrameri.cm.lapaz, 
        grid_path=None):
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
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid()

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
    a = int(np.log10(np.nanmax(MKE)))
    b = a-2
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], MKE[-1,:,:], cmap, norm, 'MKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"mke_{start_date}_{end_date}.png")
    plt.close(fig)