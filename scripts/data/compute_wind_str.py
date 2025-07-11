# LIBRAIRIES
import numpy as np
import xarray as xr
import cmocean
import glob
import os
import re
import time

store = '/lus/store/CT1/c1601279/rguillermin'
work = '/lus/work/CT1/c1601279/rguillermin'
stoens = [
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_CD', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2017_ini', 
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2018_ini_concat',
    '/lus/store/CT1/c1601279/rguillermin/run_lisa/run_swio2_stoens30_2019_ini',
    '/lus/store/CT1/c1601279/rguillermin/run_croco/run_swio2_stoens30_gls_2017',
    '/lus/store/CT1/c1601279/rguillermin/run_croco/run_swio2_stoens30_gls_2018',
    '/lus/store/CT1/c1601279/rguillermin/run_croco/run_swio2_stoens30_gls_2019'
]

grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

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
lon = g.lon_rho.isel(eta_rho = slice(0,501), xi_rho = slice(0,431))
lat = g.lat_rho.isel(eta_rho = slice(0,501), xi_rho = slice(0,431))
eta = g.eta_rho.isel(eta_rho = slice(1,502))
xi = g.xi_rho.isel(xi_rho = slice(1,432))
h = g.h
mask_rho = g.mask_rho.isel(eta_rho = slice(0,501), xi_rho = slice(0,431))
g.close()

dict = ensemble_gls_dict
for key in sorted(k for k in dict.keys() if k > '000'):
    start = time.time()
    print(f'Working on key {key}.')
    str_results = []
    
    for path in dict[key]:
        ds = xr.open_dataset(path)[['sustr', 'svstr']]
        sustr = ds.sustr.isel(eta_rho = slice(0,501)).data
        svstr = ds.svstr.isel(xi_rho = slice(0,431)).data
        sustr = xr.DataArray(
            sustr,
            name='sustr',
            dims=["time", "eta_rho", "xi_rho"],
            coords={
                "lon_rho":(["eta_rho", "xi_rho"], lon.data),
                "lat_rho":(["eta_rho", "xi_rho"], lat.data),
                "time":ds.time
            },
            attrs={
                "long_name":"averaged Kinematic u wind stress component",
                "units":"N/m2",
                "standard_name":"surface_downward_eastward_stress",
            },
        )
        svstr = xr.DataArray(
            svstr,
            name='svstr',
            dims=["time", "eta_rho", "xi_rho"],
            coords={
                "lon_rho":(["eta_rho", "xi_rho"], lon.data),
                "lat_rho":(["eta_rho", "xi_rho"], lat.data),
                "time":ds.time
            },
            attrs={
                "long_name":"averaged Kinematic v wind stress component",
                "units":"N/m2",
                "standard_name":"surface_downward_northward_stress",
            },
        )
        ds.close()
        ds = xr.merge([sustr, svstr])
        ds.attrs = {}
        ds = ds.where((mask_rho), np.nan)
        windstr = xr.DataArray(
            np.sqrt(ds.sustr**2 + ds.svstr**2),
            name='windstr',
            dims=["time", "eta_rho", "xi_rho"],
            coords={
                "lon_rho":(["eta_rho", "xi_rho"], lon.data),
                "lat_rho":(["eta_rho", "xi_rho"], lat.data),
                "time":ds.time
            },
            attrs={
                "long_name":"averaged Kinematic wind stress component",
                "units":"N/m2",
                "standard_name":"surface_downward_stress",
            },
        )
        ds = xr.merge([ds, windstr])
        str_results.append(ds)
    
    windstr = xr.concat(str_results, dim='time')
    
    path = os.path.join('/lus/work/CT1/c1601279/rguillermin/WINDSTR','run_swio2_stoens30_gls', f'{key}swiose_str.nc')
    windstr.to_netcdf(path)
    print(f'Saved {path} in {time.time() - start:.2f} s.')
