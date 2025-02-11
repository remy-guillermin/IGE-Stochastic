# utils.py
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