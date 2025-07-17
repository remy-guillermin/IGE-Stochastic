import xarray as xr
import numpy as np
import os
import glob
import sys
import time

# Path
store = '/lus/store/CT1/c1601279/lweiss/run_croco/'
scratch = '/lus/scratch/CT1/c1601279/rguillermin/CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
stoens = ['run_swio2_stoens30_2017_ini', 'run_swio2_stoens30_2017_CD', 'run_swio2_stoens30_gls_2019']
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'
NaN_output = '/lus/work/CT1/c1601279/rguillermin/NaN_CORRECTED'

ensemble_ini_files = glob.glob(os.path.join(store, 'SWIOSE_dev', stoens[0], "*swiose_avg.nc"))
print(f'{len(ensemble_ini_files)} files found for ini.')
ensemble_str_files = glob.glob(os.path.join(store, 'SWIOSE_dev', stoens[1], "*swiose_avg.nc"))
print(f'{len(ensemble_str_files)} files found for str.')
ensemble_gls_files = glob.glob(os.path.join(scratch, 'SWIOSE_dev', stoens[2], "*swiose_avg.nc"))
print(f'{len(ensemble_gls_files)} files found for gls.')

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho']]

for ensemble in stoens:
    os.makedirs(os.path.join(NaN_output, ensemble), exist_ok=True)

for f in ensemble_ini_files:
    start = time.time()
    ds = xr.open_dataset(f)[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
    ds = ds.where((g.mask_rho), np.nan)
    path = os.path.join(NaN_output, stoens[0], f.split('/')[-1])
    ds.to_netcdf(path)
    print(f"Saved {path} in {time.time() - start:.2f} s")

for f in ensemble_str_files:
    start = time.time()
    ds = xr.open_dataset(f)[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
    ds = ds.where((g.mask_rho), np.nan)
    path = os.path.join(NaN_output, stoens[1], f.split('/')[-1])
    ds.to_netcdf(path)
    print(f"Saved {path} in {time.time() - start:.2f} s")

for f in ensemble_gls_files:
    start = time.time()
    ds = xr.open_dataset(f)[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
    ds = ds.where((g.mask_rho), np.nan)
    path = os.path.join(NaN_output, stoens[2], f.split('/')[-1])
    ds.to_netcdf(path)
    print(f"Saved {path} in {time.time() - start:.2f} s")
