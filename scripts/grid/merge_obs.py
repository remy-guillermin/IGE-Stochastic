import glob
import xarray as xr
import os

obs = '/lus/work/CT1/c1601279/rguillermin/OBS'

variables = ['SST', 'SLA', 'SSS']
sources = ['OSTIA', 'SSALTO', 'MOG']

datasets = {var: glob.glob(os.path.join(obs, f'{var}*')) for var in variables}

for i, (var, path_list) in enumerate(datasets.items()):
    print(f'Opening {len(path_list)} datasets for {var}.')

    ds_list = []
    
    for path in path_list:
        ds_list.append(xr.open_dataset(path))

    ds = xr.concat(ds_list, dim='time')
    if 'depth' in ds.dims:
        ds = ds.isel(depth=0)
    
    print('Datasets concatenated.')

    path = os.path.join(obs, f'{var}_{sources[i]}')
    ds.to_netcdf(path)
    print(f'Dataset saved in {path}.')
