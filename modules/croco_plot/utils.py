"""
Module utils pour croco_plot.

Ce module contient des fonctions utilitaires pour le traitement et le chargement des données CROCO.
"""

import numpy as np
import xarray as xr

def load_grid(path):
    """
    Load the grid file into this iPython instance

    Parameters
    ----------
    path : str
        path to the 'swiose_grid.nc' file

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
    """
    g = xr.open_dataset(path)
    lon = g['lon_rho'][:, :]
    lat = g['lat_rho'][:, :]
    msk = g['mask_rho'][:, :]
    pm = g['pm'][:-1,:-1] 
    pn = g['pn'][:-1,:-1]
    msk_inv = np.where(msk == 0, msk, np.nan)
    angle = g['angle'][:, :]
    g.close()
    return lon, lat, pm, pn, msk, msk_inv, angle

def load_data(path, fields):
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
    tuple
        Tuple of loaded fields in the same order as requested.
    """
    d = xr.open_dataset(path)
    data = tuple(d[field][:, -1, :, :] for field in fields)
    d.close()
    return data

def transform_velocity(u, v, angle):
    """
    Transform the velocity components from the deformed grid to the geographic grid.

    Parameters
    ----------
    u_mean : array-like
        Surface velocity u component.
    v_mean : array-like
        Surface velocity v component.
    angle : array-like
        Grid angle values representing the grid's orientation.

    Returns
    -------
    tuple
        - u_geo: Transformed surface velocity u component
        - v_geo: Transformed surface velocity v component
    """
    u_geo = u[:-1,:].data * np.cos(angle[:-1,:-1]) - v[:,:-1].data * np.sin(angle[:-1,:-1])
    v_geo = u[:-1,:].data * np.sin(angle[:-1,:-1]) + v[:,:-1].data * np.cos(angle[:-1,:-1])
    return u_geo, v_geo

def calc_depth(s, Cs, hc, h):
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