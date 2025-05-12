import xarray as xr
import numpy as np
import glob
import os

grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'
sim = '/lus/work/CT1/c1601279/rguillermin/grid/swio_avg.nc'



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

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

ds = xr.open_dataset(sim)[['s_rho', 'Cs_rho', 'hc', 'eta_rho', 'xi_rho']]
depth_sigma = calc_depth(ds.s_rho, ds.Cs_rho, ds.hc, h)

depth_level = xr.Dataset(
    data_vars=dict(
        depth_level=(['s_rho', 'eta_rho', 'xi_rho'], depth_sigma),
    ),
    coords=dict(
        eta_rho=ds.eta_rho,
        xi_rho=ds.xi_rho,
        s_rho=ds.s_rho,
    )
)

h = h.where((h != 50), 0)
ds.close()

depth_level = depth_level.depth_level
depth_level.attrs = dict(
    units='m',
    long_name='Depth at RHO-points',
    standard_name='depth_at_rho_points',
    grid_mapping='ocean_sigma_coordinate',
)

depth_level.to_netcdf('/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc')