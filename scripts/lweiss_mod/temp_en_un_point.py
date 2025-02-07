import matplotlib.pyplot as plt
import xarray as xr

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

simu = 'run_swio_stogen_1mth_diff_100p_10jo1ud12'

depth_index = 23
lat_index = 249
lon_index = 195

ds = xr.open_dataset(data + simu + '_0001/001swiose_his.nc')  # , engine='h5netcdf')
temp = ds['temp'][:, depth_index, lat_index, lon_index]  # Température au point spécifié
time = ds['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
ds.close()

fig, ax = plt.subplots(figsize=(10, 6)) 
# Tracer la température contre le temps
ax.plot(time, temp, marker='.', linestyle='-', color='b') 

# Ajouter des labels et un titre
ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.set_title(f'Temperature at depth index {depth_index}, lat index {lat_index}, lon index {lon_index}')
 
# Afficher la figure
plt.grid(True)
plt.tight_layout()
plt.show()