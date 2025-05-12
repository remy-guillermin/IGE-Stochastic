
# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import glob
import os
import re
import time

store = '/lus/store/CT1/c1601279/rguillermin'
work = '/lus/work/CT1/c1601279/rguillermin'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
stoens = [
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_ini', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_ini_concat',
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_ini',
    '/lus/scratch/CT1/c1601279/rguillermin/CROCO/SWIOSE_dev/run_swio2_stoens30_gls_2017',
    '/lus/scratch/CT1/c1601279/rguillermin/CROCO/SWIOSE_dev/run_swio2_stoens30_gls_2018',
    '/lus/scratch/CT1/c1601279/rguillermin/CROCO/SWIOSE_dev/run_swio2_stoens30_gls_2019'
]

grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

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



ensemble_ini_files = sorted(np.concat([glob.glob(os.path.join(stoens[3], '*avg.nc')), glob.glob(os.path.join(stoens[4], '*avg.nc')), glob.glob(os.path.join(stoens[5], '*avg.nc'))]))
print(f'{len(ensemble_ini_files)} files found for INI.')
ensemble_str_files = sorted(np.concat([glob.glob(os.path.join(stoens[0], '*avg.nc')), glob.glob(os.path.join(stoens[1], '*avg.nc')), glob.glob(os.path.join(stoens[2], '*avg.nc'))]))
print(f'{len(ensemble_str_files)} files found for STR.')
ensemble_gls_files = sorted(np.concat([glob.glob(os.path.join(stoens[6], '*avg.nc')) , glob.glob(os.path.join(stoens[7], '*avg.nc')), glob.glob(os.path.join(stoens[8], '*avg.nc'))]))
print(f'{len(ensemble_gls_files)} files found for GLS.')


pattern = re.compile(r'/(\d{3})swiose_avg.nc$')
ensemble_ini_dict = {}
ensemble_str_dict = {}
ensemble_gls_dict = {}

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

for path in ensemble_gls_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_gls_dict.keys():
            ensemble_gls_dict[key] = []
        ensemble_gls_dict[key].append(path)

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


for key in sorted(k for k in ensemble_ini_dict.keys() if k > '000'):
    start = time.time()
    print(f'Working on key {key}.')
    MLD_results = []
    MLD_above_results = []
    
    for path in ensemble_ini_dict[key]:
        rho = xr.open_dataset(path)['rho']
        rho = rho.where((mask_rho), np.nan)
    
        rho_crit = rho.isel(s_rho=-2) + CRT1
        nan_mask = np.abs(np.isnan(rho_crit).mean(dim='time') - 1)
        MLDC = rho < rho_crit # Condition on density for MLD
        continuity = (MLDC != MLDC.roll(s_rho=-1)).roll(s_rho=1)
        continuity[{'s_rho': 0}] = False
        continuity_roll = continuity.roll(s_rho=-1)
    
        MLD_above = depth_level.expand_dims({'time': continuity.time}).where(continuity).max(dim='s_rho')
        
        MLD_above_results.append(MLD_above)
        
        MLD_below = depth_level.expand_dims({'time': continuity_roll.time}).where(continuity_roll).max(dim='s_rho')
        MLD_delta = MLD_below - MLD_above
        rho_above = rho.where(continuity).max(dim='s_rho')
        rho_below = rho.where(continuity_roll).max(dim='s_rho')
        
        rho_delta = rho_below - rho_above
        rho_delta = rho_delta.where(rho_delta>0)

        rho_crit_delta = rho_crit - rho_above
        rho_crit_delta = rho_crit_delta.where(rho_crit_delta>0)
        
        overlap = (~np.isnan(rho_delta) * ~np.isnan(rho_crit_delta))

        rho_delta = rho_delta.where(overlap, 0)
        rho_crit_delta = rho_crit_delta.where(overlap, 0)
        rho_pos = rho_crit_delta / rho_delta
        MLD_delta = MLD_delta.where(overlap, 0)
        MLD = MLD_above + MLD_delta * rho_pos
        
        MLD_results.append(MLD)
    
    MLD = xr.concat(MLD_results, dim='time').drop_vars('s_rho')
    MLD = MLD.rename({'depth_level': 'mld'})
    MLD['mld'].attrs = {'long_name':'interpolated mixed layer depth',
                        'units': 'meter'}
    
    MLD_above = xr.concat(MLD_above_results, dim='time')
    MLD_above = MLD_above.rename({'depth_level': 'mld_above'})
    MLD_above['mld_above'].attrs = {'long_name':'mixed layer depth',
                              'units': 'meter'}
    
    ds = xr.merge([MLD, MLD_above])
    
    path = os.path.join('/lus/work/CT1/c1601279/rguillermin/MLD','run_swio2_stoens30_ini', f'{key}swiose_mld.nc')
    ds.to_netcdf(path)
    print(f'Saved {path} in {time.time() - start:.2f} s.')
