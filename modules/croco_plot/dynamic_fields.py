"""
Module for data visualization in CROCO.

This module provides functions to plot dynamical variables (e.g., velocity, energy, vorticity) 
from CROCO simulation data or GLORYS datasets.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import cartopy.crs as ccrs
from .utils import load_grid, load_data, save_figure, plot_map

def velocity(
    data_path: str, 
    date: str, 
    figsize: tuple = (10, 8), 
    cmap = cmocean.cm.speed, 
    grid_path: str = None
):
    """
    Plot velocity data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (10, 8).
    cmap : colormap, optional
        Colormap for velocity, by default cmcrameri.cm.oslo.
    grid_path : str, optional
        Path to the grid data file, by default None.
    """
    # Load grid data
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    u = u[:,-1,:,:].sel(time=date).mean(dim='time')
    v = v[:,-1,:,:].sel(time=date).mean(dim='time')
    print("Data sliced")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    print("NaN values added")
    
    u_geo = u[:-1,:] * np.cos(angle[:-1,:-1]) - v[:,:-1] * np.sin(angle[:-1,:-1])
    v_geo = u[:-1,:] * np.sin(angle[:-1,:-1]) + v[:,:-1] * np.cos(angle[:-1,:-1])
    print("Velocity transformed")
    
    velocity = np.sqrt(u_geo**2 + v_geo**2)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Velocity SWIO {date}")
    a = -1.5
    b = 0.2
    norm = mpl.colors.LogNorm(10**a, 10**b)

    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity[:,:], cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    plt.tight_layout()
    save_figure(fig, f"velocity_{date}.png")
    plt.close(fig)

def vorticity(
    data_path: str, 
    date: str,
    figsize: tuple = (10, 8), 
    cmap = cmcrameri.cm.vik, 
    grid_path: str = None
):
    """
    Plot vorticity data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (10, 8).
    cmap : colormap, optional
        Colormap for vorticity, by default cmcrameri.cm.vik.
    grid_path : str, optional
        Path to the grid data file, by default None.
    """
    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid()

    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    u = u[:,-1,:,:].sel(time=date).mean(dim='time')
    v = v[:,-1,:,:].sel(time=date).mean(dim='time')
    print("Data sliced")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    print("NaN values added")
    
    u_geo = u[:-1,:] * np.cos(angle[:-1,:-1]) - v[:,:-1] * np.sin(angle[:-1,:-1])
    v_geo = u[:-1,:] * np.sin(angle[:-1,:-1]) + v[:,:-1] * np.cos(angle[:-1,:-1])
    print("Velocity transformed")
    
    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=1) * pm[:-1,:-1]
    du_dlat = np.gradient(u_geo, axis=0) * pn[:-1,:-1]

    # Calculate vorticity
    vorticity = dv_dlon - du_dlat

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"Vorticity SWIO {date}")
    levels = np.linspace(-0.15, 0.15, 21)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], vorticity[:,:] * 3600, cmap, norm, 'Vorticity [$h^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"vorticity_{date}.png")
    plt.close(fig)
    
def eke(
    data_path: str, 
    date: str, 
    figsize: tuple = (10, 8), 
    cmap = cmcrameri.cm.devon, 
    grid_path: str = None
):
    """
    Plot EKE data on a map for a specific date.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (10, 8).
    cmap : colormap, optional
        Colormap for EKE, by default cmcrameri.cm.lapaz.
    grid_path : str, optional
        Path to the grid data file, by default None.
    """
    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, h = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, h = load_grid()

    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
        
    u = u[:,-1,:,:]
    v = v[:,-1,:,:]
    w = w[:,-1,:,:]
    
    # Moyenne annuelle
    u_yr = u.mean(dim='time')
    v_yr = v.mean(dim='time')
    w_yr = w.mean(dim='time')
    print("Yearly mean calculated")
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan)
    v = v.where((v != fill_value), np.nan)
    w = w.where((w != fill_value), np.nan)
    print("NaN values added")
    
    # Vitesse turbulente
    ut = - (u_yr - u.sel(time=date)).data[:,:,0]
    vt = - (v_yr - v.sel(time=date)).data[:,:,0]
    wt = - (w_yr - w.sel(time=date)).data[:,:,0]
    print("Turbulent velocity calculated")
    
    ut_geo = ut[:-1,:] * np.cos(angle[:-1,:-1]) - vt[:,:-1] * np.sin(angle[:-1,:-1])
    vt_geo = ut[:-1,:] * np.sin(angle[:-1,:-1]) + vt[:,:-1] * np.cos(angle[:-1,:-1])
    wt_geo = wt[:-1,:-1]
    print("Turbulent velocity transformed")
    
    EKE = 1 / 2 * (ut_geo ** 2 + vt_geo ** 2 + wt_geo ** 2)
    print("EKE calculated")
    
    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"EKE SWIO {date}")
    a = int(np.log10(np.nanmax(EKE)))
    b = a-2.2
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], EKE, cmap, norm, 'EKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"eke_{date}.png")
    plt.close(fig)

def mke(
    data_path: str, 
    figsize: tuple = (10, 8), 
    cmap = cmcrameri.cm.devon, 
    grid_path: str = None
):
    """
    Plot MKE data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    figsize : tuple, optional
        Size of the figure, by default (10, 8).
    cmap : colormap, optional
        Colormap for MKE, by default cmcrameri.cm.lapaz.
    grid_path : str, optional
        Path to the grid data file, by default None.
    """
    # Load grid data
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid()

    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
    
    u = u[:,-1,:,:]
    v = v[:,-1,:,:]
    w = w[:,-1,:,:]
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).mean(dim='time').data
    v = v.where((v != fill_value), np.nan).mean(dim='time').data
    w = w.where((w != fill_value), np.nan).mean(dim='time').data
    print("NaN values added")
    
    u_geo = u[:-1,:] * np.cos(angle[:-1,:-1]) - v[:,:-1] * np.sin(angle[:-1,:-1])
    v_geo = u[:-1,:] * np.sin(angle[:-1,:-1]) + v[:,:-1] * np.cos(angle[:-1,:-1])
    w_geo = w[:-1,:-1]
    print("Velocity transformed")
    
    # Calculate EKE
    MKE = 1 / 2 * (u_geo ** 2 + v_geo ** 2 + w_geo ** 2)
    print("MKE calculated")

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})

    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    ax.set_title(f"MKE SWIO")
    a = int(np.log10(np.nanmax(MKE)))
    b = a-2.5
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], MKE, cmap, norm, 'MKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"mke.png")
    plt.close(fig)
    
def all(
    data_path: str, 
    date: str, 
    figsize: tuple = (30, 8), 
    vel_cmap = cmocean.cm.speed, 
    vort_cmap = cmcrameri.cm.vik,
    ke_cmap = cmcrameri.cm.devon, 
    grid_path: str = None
):

    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid()
        
    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
        
    u = u[:,-1,:,:]
    v = v[:,-1,:,:]
    w = w[:,-1,:,:]  
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan)
    v = v.where((v != fill_value), np.nan)
    w = w.where((w != fill_value), np.nan)
    print("NaN values added")    
      
    # Moyenne annuelle
    u_yr = u.mean(dim='time')
    v_yr = v.mean(dim='time')
    w_yr = w.mean(dim='time')
    print("Time mean calculated")
    
    u = u[:,:,:].sel(time=date).mean(dim='time')
    v = v[:,:,:].sel(time=date).mean(dim='time')
    w = w[:,:,:].sel(time=date).mean(dim='time')
    print("Data sliced")
    
    # Vitesse turbulente
    ut = (u - u_yr).data[:,:]
    vt = (v - v_yr).data[:,:]
    wt = (w - w_yr).data[:,:]
    print("Turbulent velocity calculated")
    
    u = u.data
    v = v.data
    w = w.data
        
    u_yr = u_yr.data
    v_yr = v_yr.data
    w_yr = w_yr.data
    
    angle = angle.data
    
    u_geo = u[:-1,:] * np.cos(angle[:-1,:-1]) - v[:,:-1] * np.sin(angle[:-1,:-1])
    v_geo = u[:-1,:] * np.sin(angle[:-1,:-1]) + v[:,:-1] * np.cos(angle[:-1,:-1])
    w_geo = w[:-1,:-1]
    
    ut_geo = ut[:-1,:] * np.cos(angle[:-1,:-1]) - vt[:,:-1] * np.sin(angle[:-1,:-1])
    vt_geo = ut[:-1,:] * np.sin(angle[:-1,:-1]) + vt[:,:-1] * np.cos(angle[:-1,:-1])
    wt_geo = wt[:-1,:-1]
    
    u_yr_geo = u_yr[:-1,:] * np.cos(angle[:-1,:-1]) - v_yr[:,:-1] * np.sin(angle[:-1,:-1])
    v_yr_geo = u_yr[:-1,:] * np.sin(angle[:-1,:-1]) + v_yr[:,:-1] * np.cos(angle[:-1,:-1])
    w_yr_geo = w_yr[:-1,:-1]
    print("Velocity transformed")
    
    velocity = np.sqrt(u_geo ** 2 + v_geo ** 2)
    velocity_t = np.sqrt(ut_geo ** 2 + vt_geo ** 2)
    velocity_yr = np.sqrt(u_yr_geo ** 2 + v_yr_geo ** 2)
    print('Velocity calculated')
    
    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=1) * pm[:-1,:-1]
    du_dlat = np.gradient(u_geo, axis=0) * pn[:-1,:-1]

    # Calculate vorticity
    vorticity = dv_dlon - du_dlat

    dvt_dlon = np.gradient(vt_geo, axis=1) * pm[:-1,:-1]
    dut_dlat = np.gradient(ut_geo, axis=0) * pn[:-1,:-1]

    vorticity_t = dvt_dlon - dut_dlat

    dv_yr_dlon = np.gradient(v_yr_geo, axis=1) * pm[:-1,:-1]
    du_yr_dlat = np.gradient(u_yr_geo, axis=0) * pn[:-1,:-1]

    vorticity_yr = dv_yr_dlon - du_yr_dlat
    print('Vorticity calculated')
    
    KE = 1 / 2 * ( u_geo ** 2 + v_geo ** 2 + w_geo ** 2)
    MKE = 1 / 2 * ( u_yr_geo ** 2 + v_yr_geo ** 2 + w_yr_geo ** 2)    
    EKE = 1 / 2 * ( ut_geo ** 2 + vt_geo ** 2 + wt_geo ** 2)
    print('Energy calculated')
    
    # Define common gridline styles
    gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

    # Plotting
    fig, axes = plt.subplots(1,3,figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    
    fig.suptitle(f"SWIO {date}")
    
    ax = axes[0]
    ax.set_title(f"Velocity")
    
    a = -1.5
    b = 0.2
    norm = mpl.colors.LogNorm(10**a, 10**b)
    
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity[:,:], vel_cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[1]
    ax.set_title(f"Vorticity")
    
    levels = np.linspace(-0.15, 0.15, 21)
    norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], vorticity[:,:] * 3600, vort_cmap, norm, 'Vorticity [$h^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[2]
    ax.set_title(f"Energy")
    
    a = int(np.log10(np.nanmax(KE)))
    b = a-2.5
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], KE, ke_cmap, norm, 'KE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"all_{date}.png")
    plt.close(fig)
    
    fig, axes = plt.subplots(1,3,figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    
    fig.suptitle(f"SWIO")
    
    ax = axes[0]
    ax.set_title(f"Velocity")
    
    a = -1.5
    b = 0.2
    norm = mpl.colors.LogNorm(10**a, 10**b)
    
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity_yr[:,:], vel_cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[1]
    ax.set_title(f"Vorticity")
    
    levels = np.linspace(-0.15, 0.15, 21)
    norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], vorticity_yr[:,:] * 3600, vort_cmap, norm, 'Vorticity [$h^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[2]
    ax.set_title(f"Energy")
    
    a = int(np.log10(np.nanmax(MKE)))
    b = a-2.5
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], MKE, ke_cmap, norm, 'MKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"all_mean.png")
    plt.close(fig)
    
    fig, axes = plt.subplots(1,3,figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    
    fig.suptitle(f"SWIO turbulent {date}")
    
    ax = axes[0]
    ax.set_title(f"Velocity")
    
    a = -1.5
    b = 0.2
    norm = mpl.colors.LogNorm(10**a, 10**b)
    
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity_t[:,:], vel_cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[1]
    ax.set_title(f"Vorticity")
    
    levels = np.linspace(-0.15, 0.15, 21)
    norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], vorticity_t[:,:] * 3600, vort_cmap, norm, 'Vorticity [$h^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)
    
    ax = axes[2]
    ax.set_title(f"Energy")
    
    a = int(np.log10(np.nanmax(EKE)))
    b = a-2.5
    norm = mpl.colors.LogNorm(vmin=10**b, vmax=10**a)
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], EKE, ke_cmap, norm, 'EKE [$m^2.s^{-2}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    save_figure(fig, f"all_turbulent_{date}.png")
    plt.close(fig)
