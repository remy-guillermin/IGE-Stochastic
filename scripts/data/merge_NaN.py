
import xarray as xr
import numpy as np
import glob
import os

work = '/lus/work/CT1/c1601279/rguillermin/'
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'
NaN_input = '/lus/work/CT1/c1601279/rguillermin/NaN_CORRECTED'
NaN_output = '/lus/work/CT1/c1601279/rguillermin/NaN_MERGED'
stoens = ['run_swio2_stoens30_ini', 'run_swio2_stoens30_str', 'run_swio2_stoens30_gls']

keys = np.array([f'{i:03}' for i in range(1,31)])

os.makedirs(os.path.join(NaN_output), exist_ok=True)
for ensemble in stoens:
    os.makedirs(os.path.join(NaN_output, ensemble), exist_ok=True)

for ensemble in stoens:
    print(f'├──Working on {ensemble}.')
    for key in keys:   
        print(f'│   ├──Working on key {key}.')
        files = glob.glob(os.path.join(NaN_input, f"*{ensemble.split('_')[-1]}", f'{key}*.nc'))
        print(f'│   │   ├──{len(files)} files found.')
        ds = xr.concat([xr.open_dataset(file) for file in files], dim='time')
        new_attrs = {
            'avg_file': ds.attrs['avg_file'],
            'grd_file': ds.attrs['grd_file'],
            'ini_file': ds.attrs['ini_file'],
            'frc_file': ds.attrs['frc_file'],    
        }
        ds.attrs = new_attrs
        path = os.path.join(NaN_output, ensemble, f'{key}swiose_avg.nc')
        ds.to_netcdf(path)
        print(f'│   │   └──File saved as {path}.')
        ds.close()
    print('│   └──Done.')
print('└──Finished.')
