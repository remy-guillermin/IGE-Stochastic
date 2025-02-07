import matplotlib.pyplot as plt
import xarray as xr

scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swiose_0002'

# Définir les indices du point spécifique
depth_index = 23
lat_index = 249
lon_index = 195

# Charger les données
ds = xr.open_dataset(scratch + simu + '/swiose_his.nc')  # , engine='h5netcdf')
temp = ds['temp'][:, depth_index, lat_index, lon_index]  # Température au point spécifié
time = ds['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
ds.close()

# Convertir les temps en format compréhensible si nécessaire (par exemple, en utilisant pandas)
# time = pd.to_datetime(time, unit='s')  # Ajustez le format de temps si nécessaire

# Créer la figure
fig, ax = plt.subplots(figsize=(10, 6))

# Tracer la température contre le temps
ax.plot(time, temp, marker='o', linestyle='-', color='b')

# Ajouter des labels et un titre
ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.set_title(f'Temperature at depth index {depth_index}, lat index {lat_index}, lon index {lon_index}')

# Afficher la figure
plt.grid(True)
plt.tight_layout()
plt.show()

#plt.savefig(path_fig + '/temp_' + str(depth_index) + '_' + str(lat_index) + '_' + str(lon_index) + '_' + simu + '.png', dpi=300, bbox_inches='tight')