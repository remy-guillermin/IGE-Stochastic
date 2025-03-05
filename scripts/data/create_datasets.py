import numpy as np
import glob
import os
import xarray as xr
import pandas as pd

#path = '/Volumes/Elements'
path = '/home/guilremy/IGE-Stochastic/data/RAW'
#output_dir = '/Users/remyguillermin/Programmation/Stage/IGE-Stochastic/data'
output_dir = '/home/guilremy/IGE-Stochastic/data/PROCESSED'
sss_path = os.path.join(path, 'SSS_CCI') # On lisa's HDD
sst_path = os.path.join(path, 'OSTIA_SST')
sla_path = os.path.join(path, 'CMEMS_SLA') # os.path.join(path, 'SLA_CMEMS')


boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)]
names = ['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']

date = '*'

sss_files = glob.glob(sss_path + f'/*{date}*.nc')
sss_files.sort()
sst_files = glob.glob(sss_path + f'/*{date}*.nc')
sst_files.sort()
sla_files = glob.glob(sss_path + f'/*{date}*.nc')
sla_files.sort()

sss_values = {name: [] for name in names}
sss_errors = {name: [] for name in names}
sst_values = {name: [] for name in names}
sst_errors = {name: [] for name in names}
sla_values = {name: [] for name in names}
sla_errors = {name: [] for name in names}
time = np.arange(np.datetime64('2017-01-01'), np.datetime64('2022-01-01'), dtype='datetime64[D]')

for file, i in zip(sss_files, len(sss_files)):
    print(f'Opening {i+1}/{len(sss_files)}')
    f = xr.open_dataset(file)
    lon = f.lon
    lat = f.lat
    sss = f.sss
    sss_err = f.sss_random_error
    f.close
    
    for name, box in zip(names, boxes):
        sss_box = sss.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        sss_err_box = sss_err.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        sss_values[name].append(np.nanmean(sss_box.values))
        sss_errors[name].append(np.nanmean(sss_err_box.values))
        
for file, i in zip(sst_files, len(sst_files)):
    print(f'Opening {i+1}/{len(sss_files)}')
    f = xr.open_dataset(file)
    lon = f.lon
    lat = f.lat
    sst = f.analysed_sst
    sst_err = f.analysis_error
    f.close
    
    for name, box in zip(names, boxes):
        sst_box = sst.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        sst_err_box = sst_err.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        sst_values[name].append(np.nanmean(sst_box.values))
        sst_errors[name].append(np.nanmean(sst_err_box.values))
        
for file, i in zip(sla_files, len(sla_files)):
    print(f'Opening {i+1}/{len(sss_files)}')
    f = xr.open_dataset(file)
    lon = f.longitude
    lat = f.latitude
    sla = f.sla
    sla_err = f.err_sla
    f.close
    
    for name, box in zip(names, boxes):
        sla_box = sla.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
        sla_err_box = sla_err.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
        sla_values[name].append(np.nanmean(sla_box.values))
        sla_errors[name].append(np.nanmean(sla_err_box.values))
        
sta_values = {name: np.array(sst_values[name]) - np.mean(sst_values[name]) for name in names}
ssa_values = {name: np.array(sss_values[name]) - np.mean(sss_values[name]) for name in names}


for name in names:
    ds = xr.Dataset({
            "sss": (["time"], sss_values[name], {
                "units": "PSU",
                "long_name": "Sea Surface Salinity",
                "standard_name": "sea_surface_salinity",
                "description": "Box averaged Sea Surface Salinity",
                "source": "CCI",
            }),
            "sss_err": (["time"], sss_errors[name], {
                "units": "PSU",
                "long_name": "Sea Surface Salinity Error",
                "description": "Uncertainty in SSS estimation",
                "source": "CCI",
            }),
            "ssa": (["time"], ssa_values[name], {
                "units": "PSU",
                "long_name": "Sea Salinity Anomaly",
                "standard_name": "sea_salinity_anomaly",
                "description": "Box averaged Sea Salinity Anomaly",
                "source": "CCI",
            }),
            "sst": (["time"], np.array(sst_values[name]) - 273.15, {  # Convert Kelvin to Celsius
                "units": "°C",
                "long_name": "Sea Surface Temperature",
                "standard_name": "sea_surface_temperature",
                "description": "Box averaged Sea Surface Temperature",
                "source": "OSTIA",
            }),
            "sst_err": (["time"], sst_errors[name], {
                "units": "°C",
                "long_name": "Sea Surface Temperature Error",
                "description": "Uncertainty in SST estimation",
                "source": "OSTIA",
            }), 
            "sta": (["time"], sta_values[name], {
                "units": "°C",
                "long_name": "Sea Temperature Anomaly",
                "standard_name": "sea_temperature_anomaly",
                "description": "Box averaged Sea Temperature Anomaly",
                "source": "OSTIA",
            }),
            "sla": (["time"], sla_values[name], { 
                "units": "m",
                "long_name": "Sea Level Anomaly",
                "standard_name": "sea_level_anomaly",
                "description": "Box averaged Sea Level Anomaly",
                "source": "CMEMS",
            }),
            "sla_err": (["time"], sla_errors[name], {
                "units": "m",
                "long_name": "Sea Level Anomaly Error",
                "description": "Uncertainty in SLA estimation",
                "source": "CMEMS",
            }),
        },
        coords={"time": time},
        attrs={
            "description": "Box averaged Sea Surface components (SSS and SST)",
            "creation_date": str(pd.Timestamp.now()),
        }
    )

    file_path = os.path.join(output_dir, f"box_data_{name}.nc")
    ds.to_netcdf(file_path)
    print(f'File {file_path} saved')