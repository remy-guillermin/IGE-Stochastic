# LIBRAIRIES
import numpy as np
import xarray as xr
import glob
import os
import re


store = '/lus/store/CT1/c1601279/rguillermin'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
stoens = [
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_ini', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_ini_concat',
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_ini'
]

grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

lat_index = [430, 280]
SWIO = (25, 69, -36, 7)

CRT1 = 0.03 # kg/m3



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



ensemble_ini_files = np.concat([glob.glob(os.path.join(stoens[3], '*avg.nc')), glob.glob(os.path.join(stoens[4], '*avg.nc')), glob.glob(os.path.join(stoens[5], '*avg.nc'))])
print(f'{len(ensemble_ini_files)} files found for INI.')
ensemble_str_files = np.concat([glob.glob(os.path.join(stoens[0], '*avg.nc')), glob.glob(os.path.join(stoens[1], '*avg.nc')), glob.glob(os.path.join(stoens[2], '*avg.nc'))])
print(f'{len(ensemble_str_files)} files found for STR.')



pattern = re.compile(r'/(\d{3})swiose_avg.nc$')
ensemble_ini_dict = {}
ensemble_str_dict = {}

for path in ensemble_ini_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_ini_dict.keys():
            ensemble_ini_dict[key] = []
        ensemble_ini_dict[key].append(path)

for path in ensemble_str_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_str_dict.keys():
            ensemble_str_dict[key] = []
        ensemble_str_dict[key].append(path)



g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
mask_rho = g.mask_rho
g.close()

ds = xr.open_dataset(ensemble_ini_files[0])[['s_rho', 'Cs_rho', 'hc', 'rho']]
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
ds.close()

h = h.where((h != 50), 0)



for key in sorted(k for k in ensemble_str_dict.keys() if k > '016'):
    print(f'Working on key {key}.')
    MLD_results = []
    for path in ensemble_str_dict[key]:
        ds = xr.open_dataset(path)['rho']
        ds = ds.where((mask_rho), np.nan)
        
        MLDC = ds - ds.isel(s_rho=-2) < CRT1 # Condition on density for MLD
        mask = - MLDC.sum(dim='s_rho') 
        MLD = depth_level.isel(s_rho=mask)
        MLD = MLD.where((g.mask_rho), np.nan)

        MLD_results.append(MLD)
        
    MLD = xr.concat(MLD_results, dim='time').drop_vars('s_rho')
    MLD = MLD.rename({'depth_level': 'mld'})
    MLD['mld'].attrs = {'long_name':'mixed layer depth',
                        'units': 'meter'}
    path = os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD','run_swio2_stoens30_CD', f'{key}swiose_mld.nc')
    MLD.to_netcdf(path)
    print(f'Saved {path}.')



# for key in sorted(k for k in ensemble_ini_dict.keys() if k > '029'):
#     print(f'Working on key {key}.')
#     MLD_results = []
#     for path in ensemble_ini_dict[key]:
#         ds = xr.open_dataset(path)['rho']
#         ds = ds.where((mask_rho), np.nan)
        
#         MLDC = ds - ds.isel(s_rho=-2) < CRT1 # Condition on density for MLD
#         mask = - MLDC.sum(dim='s_rho') 
#         MLD = depth_level.isel(s_rho=mask)
#         MLD = MLD.where((g.mask_rho), np.nan)

#         MLD_results.append(MLD)
        
#     MLD = xr.concat(MLD_results, dim='time').drop_vars('s_rho')
#     MLD = MLD.rename({'depth_level': 'mld'})
#     MLD['mld'].attrs = {'long_name':'mixed layer depth',
#                         'units': 'meter'}
#     path = os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD','run_swio2_stoens30_ini', f'{key}swiose_mld.nc')
#     MLD.to_netcdf(path)
#     print(f'Saved {path}.')






