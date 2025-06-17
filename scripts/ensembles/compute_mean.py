import xarray as xr
import glob
import os

files = glob.glob(os.path.join('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_gls', '0*.nc'))
files.sort()
print('Found %d files' % len(files))

ds_concat = xr.concat([xr.open_dataset(file) for file in files], dim='member')
print('Concatenation done')

mean = ds_concat.mean(dim='member')
mean.to_netcdf('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_gls/mean_swiose_mld.nc')
std = ds_concat.std(dim='member')
std.to_netcdf('/lus/work/CT1/c1601279/rguillermin/NaN_MERGED/run_swio2_stoens30_gls/std_swiose_mld.nc')
