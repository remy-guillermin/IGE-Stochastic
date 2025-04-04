import xarray as xr
import numpy as np
import os
import glob
import sys
import time

# Path
store = '/lus/store/CT1/c1601279/lweiss/run_croco/'
scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE_dev/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
stoens = ['run_swio2_stoens30_2017_ini', 'run_swio2_stoens30_2017_CD']
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
NaN_output = '/lus/work/CT1/c1601279/rguillermin/NaN_CORRECTED'

ensemble_ini_files = [os.path.join(scratch, stoens[0], f"{i:03d}swiose_avg.nc") for i in range(1, 31)]
ensemble_str_files = [os.path.join(store, 'SWIOSE_dev', stoens[1], f"{i:03d}swiose_avg.nc") for i in range(1, 31)]

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
