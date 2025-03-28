# Libs
import xarray as xr
import xesmf as xe
from glob import glob
import os
import time

# Need to install xesmf with conda: conda install -c conda-forge xesmf (because of the dependency of esmpy that cannot be installed with pip) for LOCAL
# Need to manually install esmpy then install xesmf with pip for CLUSTER (see README.md)

glob_start = time.time()

# LOCAL
# grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
# data_path = '/home/guilremy/IGE-Stochastic/Data/data/RAW'

# CLUSTER
grid_path = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
data_path = '/lus/work/CT1/c1601279/rguillermin/OBS/'

grid = xr.open_dataset(grid_path)
lon = grid.lon_rho.values
lat = grid.lat_rho.values
grid.close()

print('Grid loaded')

grid_out = xr.Dataset(
    {
        "lat" : (["y", "x"], lat),
        "lon" : (["y", "x"], lon),
    }
)

# data = glob(os.path.join(data_path, 'SSS_CMEMS', '*.nc'))     # SSS
# data = glob(os.path.join(data_path, 'SLA_CMEMS', '*.nc'))     # SLA
data = glob(os.path.join(data_path, 'SST_OSTIA', '*.nc'))     # SST
data.sort()

path = data[0]

#ds = xr.open_dataset(path)
ds = xr.open_dataset(path)[['analysed_sst', 'analysis_error', 'sea_ice_fraction']]

path = data[0]

ds = xr.open_dataset(path)

regridder_start = time.time()
regridder = xe.Regridder(ds, grid_out, 'bilinear')

print(f'Regridder created in {time.time() - regridder_start:.2f} s')

datasets = []

for path in data[:1000]:
    file_start = time.time()
    ds = xr.open_dataset(path)[['analysed_sst', 'analysis_error', 'sea_ice_fraction']]
    ds_regridded = regridder(ds)
    ds.close()
    # datasets.append(ds_regridded[['sos', 'sos_error']])                 # SSS
    # datasets.append(ds_regridded[['sla', 'err_sla']])                   # SLA  
    datasets.append(ds_regridded[['analysed_sst', 'analysis_error']])   # SST
    print(f'{path} regridded in {time.time() - file_start:.2f} s')

ds_combined = xr.concat(datasets, dim='time')

# output_path = '/home/guilremy/IGE-Stochastic/Data/data/REGRIDDED/SST_OSTIA_part1.nc'      # On local
output_path = '/lus/work/CT1/c1601279/rguillermin/OBS/SST_OSTIA_part1.nc'                 # On cluster
ds_combined.to_netcdf(output_path)
print(f'Saved in {output_path} in {time.time() - glob_start:.2f} s')
