"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage des données CROCO spécifiquement dans la région du courant des Aiguilles.
"""
import os
from pathlib import Path
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from .utils import load_grid, load_data, save_figure, plot_zoom, plot_map

def zoom_velocity(
    data_path, 
    date,
    figsize=(10, 8), 
    cmap=cmocean.cm.speed, 
    grid_path=None,
    isFilm=False
):
    """
    Plot velocity data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for velocity, by default cmcrameri.cm.oslo
    """
    if isFilm:
        output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
        folder_path = Path(output_dir) / Path(os.path.splitext(data_path)[0])
        # Create the folder if it doesn't exist
        folder_path.mkdir(parents=True, exist_ok=True)
        print(f"Folder '{folder_path}' is ready.")

    # Load grid data
    if grid_path is not None:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, msk, msk_inv, angle, _ = load_grid()
        
    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    if date not in u['time'].dt.strftime('%Y-%m-%d').values:
        raise ValueError(f"Date {date} not found in the time dimension of u.")
    if date not in v['time'].dt.strftime('%Y-%m-%d').values:
        raise ValueError(f"Date {date} not found in the time dimension of v.")
    
    u = u[:,:,:,:].sel(time=date).mean(dim='time')
    v = v[:,:,:,:].sel(time=date).mean(dim='time')
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

    ax.set_title(f"Velocity SWIO {date}", size=9)
    a = -1.2
    b = 0
    norm = mpl.colors.LogNorm(10**a, 10**b)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    plot_zoom(ax, lon[:-1,:-1], lat[:-1,:-1], velocity[-1,:,:], cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style)

    plt.tight_layout()
    if isFilm:
        save_figure(fig, f"{os.path.splitext(data_path)[0]}/agulhas_zoom_velocity_{date}.png")
    else:
        save_figure(fig, f"agulhas_zoom_velocity_{date}.png")
    plt.close(fig)
    
def velocity(
    data_path, 
    date,
    figsize=(10, 8), 
    cmap=cmocean.cm.speed, 
    grid_path=None,
    isFilm=False
):
    """
    Plot velocity data on a map for a specific date range.

    Parameters
    ----------
    data_path : str
        Path to the simulation data file.
    start_date : str
        Date for the data slice in 'YYYY-MM-DD' format.
    figsize : tuple, optional
        Size of the figure, by default (8, 8)
    cmap : colormap, optional
        Colormap for velocity, by default cmcrameri.cm.oslo
    """
    if isFilm:
        output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
        folder_path = Path(output_dir) / Path(os.path.splitext(data_path)[0])
        # Create the folder if it doesn't exist
        folder_path.mkdir(parents=True, exist_ok=True)

        print(f"Folder '{folder_path}' is ready.")

    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid()
        
    # Load simulation data
    u, v = load_data(data_path, ('u', 'v'))
    
    if date not in u['time'].dt.strftime('%Y-%m-%d').values:
        raise ValueError(f"Date {date} not found in the time dimension of u.")
    if date not in v['time'].dt.strftime('%Y-%m-%d').values:
        raise ValueError(f"Date {date} not found in the time dimension of v.")
    
    u = u[:,:,:,:].sel(time=date).mean(dim='time')
    v = v[:,:,:,:].sel(time=date).mean(dim='time')
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

    ax.set_title(f"Velocity SWIO {date}", size=9)
    #levels = np.linspace(0, 2.5, 26)
    a = -1.2
    b = 0
    norm = mpl.colors.LogNorm(10**a, 10**b)
    
    zoom = np.min(lon.data) + 1 , 50, np.min(lat.data) + 1, - 10.0
    
    plot_map(ax, lon[:-1,:-1], lat[:-1,:-1], velocity[-1,:,:], cmap, norm, 'Velocity [$m.s^{-1}$]', msk[:-1,:-1], msk_inv[:-1,:-1], gridline_style, bounds=zoom)

    plt.tight_layout()
    if isFilm:
        save_figure(fig, f"agulhas_velocity_{date}.png", isFilm=isFilm, filmDir=os.path.splitext(data_path)[0])
    else:
        save_figure(fig, f"agulhas_velocity_{date}.png")
    plt.close(fig)