"""
Utility module for croco_plot.

This module contains utility functions for processing and loading CROCO data.
"""

import numpy as np
import xarray as xr
import os
import subprocess
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
import cartopy.crs as ccrs
from typing import Optional, Tuple, Union, Dict, List
import ffmpeg
import glob

R = 6371000.0  

def load_grid(
    path: Optional[str] = None
) -> Tuple[
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
    Load the grid file and extract relevant grid parameters.

    Parameters
    ----------
    path : str, optional
        Path to the grid file. If None, a default path is used.

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
        path = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
    g = xr.open_dataset(path)
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

def load_data(
    path: str, 
    fields: Tuple[str, ...]
) -> Union[Tuple[xr.DataArray, ...], xr.DataArray]:
    """
    Load specified fields from a simulation data file.

    Parameters
    ----------
    path : str
        Path to the simulation data file.
    fields : tuple of str
        Tuple of field names to load (e.g., ('u', 'v', 'temp', 'salt')).

    Returns
    -------
    Union[Tuple[xr.DataArray, ...], xr.DataArray]
        - If multiple fields are requested, returns a tuple of DataArray objects in the same order as requested.
        - If only one field is requested, returns a single DataArray object.
    """
    d = xr.open_dataset(path)
    if len(fields) == 1:
        data = d[fields[0]]  # Return DataArray instead of values
    else:
        data = tuple(d[field] for field in fields)
    d.close()
    print(f"Data in {path} loaded.")
    return data

def load_GLORYS(
    path: str,
) -> Tuple[xr.DataArray, ...]:
    """
    Load fields from a GLORYS data file.

    Parameters
    ----------
    path : str
        Path to the GLORYS data file.

    Returns
    -------
    tuple
        - salt: Salinity values.
        - temp: Temperature values.
        - u: Zonal velocity values.
        - v: Meridional velocity values.
        - zeta: Sea surface height values.
        - Lon: Longitude grid values.
        - Lat: Latitude grid values.
        - msk: Mask array of valid grid points.
        - msk_inv: Inverse mask array with invalid points set to NaN.
    """
    d = xr.open_dataset(path)
    salt = d['so']
    temp = d['thetao']
    zeta = d['zos']
    u = d['uo']
    v = d['vo']
    lon = d['longitude']
    lat = d['latitude']
    d.close()
    msk = ~np.isnan(salt[0, 0, :, :]).data
    msk_inv = np.where(msk == 0, msk, np.nan)
    Lon, Lat = np.meshgrid(lon, lat)
    print(f"GLORYS data in {path} loaded.")
    return salt, temp, zeta, u, v, Lon, Lat, msk, msk_inv

def calc_depth(
    s: xr.DataArray, 
    Cs: xr.DataArray, 
    hc: float, 
    h: xr.DataArray
) -> xr.DataArray:
    """
    Compute depth using the S-coordinate transformation.

    Parameters
    ----------
    s : xr.DataArray
        S-coordinate at RHO-points, typically ranging from -1 (surface) to 0 (bottom).
    Cs : xr.DataArray
        S-coordinate stretching curves at W-points, defining the vertical stretching.
    hc : float
        Critical depth parameter (in meters), influencing vertical terrain-following transformation.
    h : xr.DataArray
        Bathymetric depth at RHO-points (in meters), representing the seafloor depth.

    Returns
    -------
    xr.DataArray
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

def haversine(
    lat1, 
    lon1, 
    lat2, 
    lon2
):
    """
    Compute the Haversine distance between two points on Earth.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : float
        Latitude and longitude of the two points in degrees.

    Returns
    -------
    float
        Distance between the two points in meters.
    """
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c 

def compute_dx_dy(
    lat, 
    lon
):
    """
    Compute effective grid spacing (dx_eff and dy_eff) for a 2D grid.

    Parameters
    ----------
    lat : np.ndarray
        2D array of latitudes.
    lon : np.ndarray
        2D array of longitudes.

    Returns
    -------
    tuple
        - dx_eff: Effective zonal spacing in meters.
        - dy_eff: Effective meridional spacing in meters.
    """
    dy_eff = haversine(lat[:-1, :], lon[:-1, :], lat[1:, :], lon[1:, :])
    dx_eff = haversine(lat[:, :-1], lon[:, :-1], lat[:, 1:], lon[:, 1:])

    return dx_eff, dy_eff

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
    levels: Optional[np.ndarray] = None,
    bounds: Optional[Tuple[float, float, float, float]] = None,
    interactive: Optional[bool] = False
) -> None:
    """
    Plot data on a map with optional gridlines and bounds.

    Parameters
    ----------
    ax : Axes
        The axis to plot on.
    lon : np.ndarray
        Longitudes.
    lat : np.ndarray
        Latitudes.
    data : np.ndarray
        Data to plot.
    cmap : Colormap
        Colormap to use.
    norm : Normalize
        Normalization for the colormap.
    label : str
        Label for the colorbar.
    msk : np.ndarray
        Mask for contour.
    msk_inv : np.ndarray
        Inverse mask for contourf.
    gridline_style : dict
        Style for gridlines.
    levels : np.ndarray, optional
        Contour levels.
    bounds : tuple of float, optional
        Bounds for the plot (x1, x2, y1, y2).
    interactive : bool, optional
        Whether to use an interactive backend for plotting.
    """
    if interactive:
        original_backend = matplotlib.get_backend()  # Store the original backend
        matplotlib.use('tkagg')  # Temporarily switch to interactive backend
    try:
        pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        ax.contourf(lon, lat, msk_inv, colors='lightgray')
        ax.contour(lon, lat, msk, colors='k', linewidths=0.2)
        
        if bounds is not None:
            x1, x2, y1, y2 = bounds
            ax.set_xlim(x1, x2)
            ax.set_ylim(y1, y2)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'color': 'k'}

        cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
        if levels is not None:
            ticks = levels
            cb.set_ticks(ticks)
            cb.ax.set_yticklabels(np.round(ticks, 2))
    except Exception as e:
        print(f"Error in plot_map: {e}")
    finally:
        if interactive:
            matplotlib.use(original_backend)  # Restore the original backend after processing
 
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
    levels: Optional[np.ndarray] = None,
    interactive: Optional[bool] = False
) -> None:
    """
    Plot data with a zoomed inset on a map.

    Parameters
    ----------
    ax : Axes
        The axis to plot on.
    lon : np.ndarray
        Longitudes.
    lat : np.ndarray
        Latitudes.
    data : np.ndarray
        Data to plot.
    cmap : Colormap
        Colormap to use.
    norm : Normalize
        Normalization for the colormap.
    label : str
        Label for the colorbar.
    msk : np.ndarray
        Mask for contour.
    msk_inv : np.ndarray
        Inverse mask for contourf.
    gridline_style : dict
        Style for gridlines.
    levels : np.ndarray, optional
        Contour levels.
    interactive : bool, optional
        Whether to use an interactive backend for plotting.
    """
    if interactive:
        original_backend = matplotlib.get_backend()  # Store the original backend
        matplotlib.use('tkagg')  # Temporarily switch to interactive backend
    try:
        x1, x2, y1, y2 = np.min(lon.data) + 1 , 50, np.min(lat.data) + 1, - 10.0

        pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        ax.contourf(lon, lat, msk_inv, colors='lightgray')
        ax.contour(lon, lat, msk, colors='k', linewidths=0.2)
        
        axins = ax.inset_axes([0.33, 0.03, 0.64, 0.94], projection=ccrs.PlateCarree(), anchor='NE')
        axins.set_extent(ax.get_extent(), crs=ccrs.PlateCarree())

        # Inset map using pcolormesh
        axins.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        axins.contourf(lon, lat, msk_inv, colors='lightgray')
        axins.contour(lon, lat, msk, colors='k', linewidths=0.2)
        
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'color': 'k'}
        
        ax.plot([x1, x2], [y1, y1], "k--")
        ax.plot([x2, x2], [y1, y2], "k--")
        ax.plot([x2, x1], [y2, y2], "k--")
        ax.plot([x1, x1], [y2, y1], "k--")
        
        cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
        if levels is not None:
            ticks = levels
            cb.set_ticks(ticks)
            cb.ax.set_yticklabels(np.round(ticks, 2))
    except Exception as e:
        print(f"Error in plot_zoom: {e}") 
    finally:
        matplotlib.use(original_backend) 
        
def plot_time_series(
    time_results: np.ndarray,
    results: Dict[str, np.ndarray],
    ylabel: str,
    roll: int,
    names: List[str],
    colors: List[str],
    filename: str,
    title: str,
    interactive: Optional[bool] = False
) -> None:
    """
    Plot time series data with rolling mean for multiple variables.

    Parameters
    ----------
    time_results : np.ndarray
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
    interactive : bool, optional
        Whether to use an interactive backend for plotting.
    """
    if interactive:
        original_backend = matplotlib.get_backend()  # Store the original backend
        matplotlib.use('tkagg')  # Temporarily switch to interactive backend
    try:
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
    except Exception as e:
        print(f"Error in plot_time_series: {e}") 
    finally:
        matplotlib.use(original_backend)  # Restore the original backend after processing

def prepare_film(
    func: callable, 
    data_path: str, 
    dates: List[str], 
    grid_path: Optional[str] = None, 
    **kwargs
) -> None:
    """
    Prepare a series of plots for a film by varying the date parameter.

    Parameters
    ----------
    func : callable
        The plotting function to use.
    data_path : str
        Path to the simulation data file.
    dates : list of str
        List of dates for the data slices in 'YYYY-MM-DD' format.
    grid_path : str, optional
        Path to the grid data file, by default None.
    **kwargs : dict
        Additional keyword arguments to pass to the plotting function.
    """
    original_backend = matplotlib.get_backend()  # Store the original backend
    matplotlib.use('Agg')  # Temporarily switch to non-interactive backend

    try:
        for date in dates:
            print(f"Processing date: {date}")
            func(data_path=data_path, date=date, grid_path=grid_path, isFilm=True, **kwargs)
    finally:
        matplotlib.use(original_backend)  # Restore the original backend after processing
        
def create_film(
    filmDir: str, 
    filmName: str,
    framerate: Optional[int] = 25,
) -> None:
    """
    Create a film from saved figures.

    Parameters
    ----------
    filmDir : str
        Directory containing the film frames.
    filmName : str
        Name of the film file to create, which is the same as the frame name before the date stamp.
    framerate : int, optional
        Framerate of the output video (default is 25 fps).
    """
    output_dir = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
    film_path = os.path.join(output_dir, filmDir)
    film_file = os.path.join(output_dir, filmName)
    images = glob.glob(f'{film_path}/{filmName}_*.png')
    images.sort()
    input_files = '|'.join(images)
    
    ffmpeg.input(f'concat:{input_files}', r=framerate).output(f'{film_file}.mp4', pix_fmt='yuv420p').run(overwrite_output=True)
    
    print(f"""
Film created as {film_file}.mp4.
""")

def save_figure(
    fig: plt.Figure, 
    filename: str,
    filmDir: Optional[str] = None,
    isFilm: Optional[bool] = False,
    isTransparent: Optional[bool] = True
) -> None:
    """
    Save the figure to the specified filename.

    Parameters
    ----------
    fig : plt.Figure
        The figure to save.
    filename : str
        The path to save the figure.
    filmDir : str, optional
        Directory for film frames, required if isFilm is True.
    isFilm : bool, optional
        Whether the figure is part of a film sequence, by default False.
    isTransparent : bool, optional
        Whether to save the figure with a transparent background, by default True.
    """
    home_dir = os.path.expanduser("~") # Get the home directory
    output_dir = os.path.join(home_dir, "Images")
    if isFilm and filmDir is not None:
        output_dir = os.path.join(output_dir, filmDir) 
    elif isFilm and filmDir is None:
        raise ValueError(f"Film directory (filmDir) is not specified")
    os.makedirs(output_dir, exist_ok=True)
    
    fig.savefig(os.path.join(output_dir, filename), dpi=300, transparent=isTransparent)
    if isFilm:
        print(f"Figure saved as {os.path.join(output_dir, filename)}.")
    else:
        print(f"""
Figure saved as {filename}.
To open the figure, run: cplot.utils.open_figure('{filename}')
""")

def open_figure(
    filenames: Union[str, List[str]]
) -> None:
    """
    Open saved figure(s) using a terminal command.

    Parameters
    ----------
    filenames : str or list of str
        The name of the file or list of files to open.
    """
    home_dir = os.path.expanduser("~") # Get the home directory
    output_dir = os.path.join(home_dir, "Images")
    
    if isinstance(filenames, str):
        filenames = [filenames]

    file_paths = [os.path.join(output_dir, filename) for filename in filenames]

    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist.")

    for file_path in file_paths:
        subprocess.Popen(['eog', file_path])  # Use eog for Linux
        
def get_color_from_filename(
    filename, 
    names:  Optional[List[str]] = ['Equator', 'Mascarene', 'Mayotte-Comores', 'South-Moz'], 
    colors: Optional[List[str]] = ['saddlebrown', 'teal', 'darkorchid', 'navy']
):
    """
    Get a color based on the filename.

    Parameters
    ----------
    filename : str
        The filename to check.
    names : list of str, optional
        List of region names to match.
    colors : list of str, optional
        List of colors corresponding to the region names.

    Returns
    -------
    str
        The color corresponding to the matched region, or 'black' if no match is found.
    """
    return_color = 'black'
    for name, color in zip(names, colors):
        if name in filename:
            return_color = color
    return return_color

def get_name_from_filename(
    filename, 
    names=['Equator', 'Mascarene', 'Mayotte-Comores', 'South-Moz']
):
    """
    Get a region name based on the filename.

    Parameters
    ----------
    filename : str
        The filename to check.
    names : list of str
        List of region names to match.

    Returns
    -------
    str
        The matched region name, or 'global' if no match is found.
    """
    return_name = 'global'
    for name in names:
        if name in filename:
            return_name = name
    return return_name
