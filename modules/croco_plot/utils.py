"""
Module utils pour croco_plot.

Ce module contient des fonctions utilitaires pour le traitement et le chargement des données CROCO.
"""

import numpy as np
import xarray as xr
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
import cartopy.crs as ccrs
from typing import Optional, Tuple, Union, Dict, List

def load_grid(path: Optional[str] = None, is_Velocity: bool = False) -> Tuple[
    np.ndarray,  # lon
    np.ndarray,  # lat
    np.ndarray,  # pm
    np.ndarray,  # pn
    np.ndarray,  # msk
    np.ndarray,  # msk_inv
    np.ndarray,  # angle
    np.ndarray   # h
]:
    """
    Load the grid file into this iPython instance

    Returns
    -------
    tuple
        - lon: Longitude grid values.
        - lat: Latitude grid values.
        - pm: Curvilinear coordinate metric in XI.
        - pn: Curvilinear coordinate metric in ETA.
        - msk: Mask array of valid grid points.
        - msk_inv: Inverse mask array with invalid points set to NaN.
        - angle: Grid angle values representing the grid's orientation.
        - h: Bathymetric depth values.
    """
    if path is None:
        path = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'
    g = xr.open_dataset(path)
    if is_Velocity:
        lon = g['lon_rho'][:-1, :-1]
        lat = g['lat_rho'][:-1, :-1]
        msk = g['mask_rho'][:-1, :-1]
        pm = g['pm'][:-1,:-1] 
        pn = g['pn'][:-1,:-1]
        msk_inv = np.where(msk == 0, msk, np.nan)
        h = g['h'][:-1, :-1]
    else:
        lon = g['lon_rho'][:, :]
        lat = g['lat_rho'][:, :]
        msk = g['mask_rho'][:, :]
        pm = g['pm'][:,:] 
        pn = g['pn'][:,:]
        msk_inv = np.where(msk == 0, msk, np.nan)
        h = g['h'][:, :]
    angle = g['angle'][:, :]
    g.close()
    print("Grid loaded.")
    return lon, lat, pm, pn, msk, msk_inv, angle, h

def load_data(path: str, fields: Tuple[str, ...]) -> Union[Tuple[xr.DataArray, ...], xr.DataArray]:
    """
    Load the specified fields from the simulation data file.

    Parameters
    ----------
    path : str
        Path to the simulation data file.
    fields : tuple
        Tuple of field names to load (e.g., ('u', 'v', 'temp', 'salt')).

    Returns
    -------
    tuple or DataArray
        Tuple of loaded fields in the same order as requested if multiple fields are requested,
        otherwise a single DataArray if only one field is requested.
    """
    d = xr.open_dataset(path)
    if len(fields) == 1:
        data = d[fields[0]]  # Return DataArray instead of values
    else:
        data = tuple(d[field] for field in fields)
    d.close()
    print(f"Data in {path} loaded.")
    return data

def calc_depth(s: xr.DataArray, Cs: xr.DataArray, hc: float, h: xr.DataArray) -> xr.DataArray:
    """
    Compute the depth of the mask using the S-coordinate transformation.

    Parameters
    ----------
    s : array-like
        S-coordinate at RHO-points, typically ranging from -1 (surface) to 0 (bottom).
    Cs : array-like
        S-coordinate stretching curves at W-points, defining the vertical stretching.
    hc : float
        Critical depth parameter (in meters), influencing vertical terrain-following transformation.
    h : array-like
        Bathymetric depth at RHO-points (in meters), representing the seafloor depth.

    Returns
    -------
    array-like
        Computed depth at RHO-points.
    """
    N = len(s)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h
    return depth

def plot_map(
    ax: Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    data: np.ndarray,
    cmap: Colormap,
    norm: Normalize,
    label: str,
    msk: np.ndarray,
    msk_inv: np.ndarray,
    gridline_style: Dict,
    levels: Optional[np.ndarray] = None
) -> None:
    """
    Helper function to plot data on a given axis.

    Parameters
    ----------
    ax : Axes
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
    try:
        pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
        ax.contourf(lon, lat, msk_inv, colors='lightgray')

        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}

        cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
        if levels is not None:
            ticks = levels
            cb.set_ticks(ticks)
            cb.ax.set_yticklabels(np.round(ticks, 2), fontsize=8)
    except Exception as e:
        print(f"Error in plot_map: {e}")
 
def plot_zoom(
    ax: Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    data: np.ndarray,
    cmap: Colormap,
    norm: Normalize,
    label: str,
    msk: np.ndarray,
    msk_inv: np.ndarray,
    gridline_style: Dict,
    levels: Optional[np.ndarray] = None
) -> None:
    """
    Helper function to plot data on a given axis.

    Parameters
    ----------
    ax : Axes
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
    try:
        x1, x2, y1, y2 = np.min(lon.data) + 1 , 40.0, np.min(lat.data) + 1, - 15.0

        pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
        ax.contourf(lon, lat, msk_inv, colors='lightgray')
        
        axins = ax.inset_axes([0.33, 0.03, 0.64, 0.94], projection=ccrs.PlateCarree(), anchor='NE')
        axins.set_extent(ax.get_extent(), crs=ccrs.PlateCarree())

        # Inset map using pcolormesh
        axins.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        axins.contour(lon, lat, msk, colors='k', linewidths=0.1)
        axins.contourf(lon, lat, msk_inv, colors='lightgray')
        
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}
        
        glins = axins.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        glins.top_labels = False
        glins.right_labels = False
        glins.bottom_labels = False
        glins.left_labels = False
        
        ax.plot([x1, x2], [y1, y1], "k--")
        ax.plot([x2, x2], [y1, y2], "k--")
        ax.plot([x2, x1], [y2, y2], "k--")
        ax.plot([x1, x1], [y2, y1], "k--")
        
        cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
        if levels is not None:
            ticks = levels
            cb.set_ticks(ticks)
            cb.ax.set_yticklabels(np.round(ticks, 2), fontsize=8)
    except Exception as e:
        print(f"Error in plot_map: {e}") 
        
def plot_time_series(
    time_results: np.ndarray,
    results: Dict[str, np.ndarray],
    ylabel: str,
    roll: int,
    names: List[str],
    colors: List[str],
    filename: str,
    title: str
) -> None:
    """
    Plot time series data with rolling mean for multiple variables.

    Parameters
    ----------
    time_results : array-like
        Array of time points corresponding to the results.
    results : dict
        Dictionary containing the results to plot, with keys as variable names and values as arrays of data.
    ylabel : str
        Label for the y-axis.
    roll : int
        Window size for computing the rolling mean.
    names : list of str
        List of variable names to plot.
    colors : list of str
        List of colors for each variable plot.
    filename : str
        Filename to save the plot.
    title : str
        Title of the plot.
    """
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)
    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, results[name], color=color, linestyle='--', linewidth=1)
        rolling_mean = np.convolve(results[name], np.ones(roll)/roll, mode='same')
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel(ylabel)
        ax.set_title(name)
    axes[-1].set_xlabel('Time')
    fig.suptitle(title)
    save_figure(fig, filename)
    plt.close(fig)

def save_figure(fig: plt.Figure, filename: str) -> None:
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
    fig.savefig(os.path.join(output_dir, filename), dpi=300)
    print(f"""
Figure saved as {filename}.
To open the figure, run: cplot.utils.open_figure('{filename}')
""")


def open_figure(filenames: Union[str, List[str]]) -> None:
    """
    Open the saved figure(s) using a terminal command without blocking the IPython session.

    Parameters
    ----------
    filenames : str or list of str
        The name of the file or list of files to open.
    """
    output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
    
    if isinstance(filenames, str):
        filenames = [filenames]

    file_paths = [os.path.join(output_dir, filename) for filename in filenames]

    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist.")

    for file_path in file_paths:
        subprocess.Popen(['eog', file_path])  # Use eog for Linux