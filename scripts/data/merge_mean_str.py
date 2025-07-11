# LIBRAIRIES
import numpy as np
import xarray as xr
import glob
import os

path = '/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_gls/'
files = glob.glob(os.path.join(path, '*str.nc'))
files.sort()
files

ds_concat = xr.concat([xr.open_dataset(file).drop_duplicates(dim='time', keep='last') for file in files], dim='member')
print('Concatenation done')

mean = ds_concat.mean(dim='member')
mean.to_netcdf('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_gls/mean_swiose_str.nc')
print('Mean saved')
std = ds_concat.std(dim='member')
std.to_netcdf('/lus/work/CT1/c1601279/rguillermin/WINDSTR/run_swio2_stoens30_gls/std_swiose_str.nc')
print('Std saved')