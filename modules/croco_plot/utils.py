"""
Module utils pour croco_plot.

Ce module contient des fonctions utilitaires pour le traitement et l'affichage
des données CROCO.
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
        - msk: Mask array of valid grid points.
        - msk_inv: Inverse mask array with invalid points set to NaN.
        - angle: Grid angle values representing the grid's orientation.
    """
    g = xr.open_dataset(path)
    lon = g['lon_rho'][:, :]
    lat = g['lat_rho'][:, :]
    msk = g['mask_rho'][:, :]
    msk_inv = np.where(msk == 0, msk, np.nan)
    # dx = 1 / g['pm'] * units.meter
    # dy = 1 / g['pn'] * units.meter
    angle = g['angle'][:, :]
    #sys.exit()
    g.close()
    return lon, lat, msk, msk_inv, angle

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