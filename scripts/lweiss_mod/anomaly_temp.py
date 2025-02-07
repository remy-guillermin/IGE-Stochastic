import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

simu = 'run_swio_stogen_1mth_diff_100p_10jo1ud12'

# Charger les données
ds = xr.open_dataset(data + simu + '_0001/001swiose_his.nc')  # , engine='h5netcdf')
temp = ds['temp'][:, :, :, :]  # Sélectionne le niveau de surface s_rho
time = ds['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
lon = ds['lon_rho'][:, :]
lat = ds['lat_rho'][:, :]
msk = ds['mask_rho'][:, :]
ds.close()

# Initialiser les listes pour stocker les températures maximales et leurs indices
max_temps = []
max_indices = []

# Parcourir chaque pas de temps pour trouver la température maximale et ses indices
for t in range(temp.shape[0]):
    temp_t = temp[t, :, :, :]
    masked_temp_t = np.where(msk == 1, temp_t, np.nan)  # Appliquer le masque
    max_temp = np.nanmin(masked_temp_t)
    max_temps.append(max_temp)
    
    # Trouver les indices de la température maximale
    max_idx = np.unravel_index(np.nanargmin(masked_temp_t), masked_temp_t.shape)
    max_indices.append(max_idx)

# Convertir les indices en numpy array pour faciliter l'accès
max_indices = np.array(max_indices)

# Créer la figure
fig, ax = plt.subplots(figsize=(10, 6))

# Tracer la température maximale contre le temps
ax.plot(time, max_temps, marker='o', linestyle='-')

# Annoter les indices de la grille de la température maximale
for i, txt in enumerate(max_indices):
    ax.annotate(f'{txt}', (time[i], max_temps[i]), textcoords="offset points", xytext=(0,10), ha='center', size=6)

# Ajouter des labels et un titre
ax.set_xlabel('Time')
ax.set_ylabel('Temperature (°C)')
plt.show()
#plt.savefig(path_fig + '/temp_minimum_anomaly_' + simu + '.png', dpi=300, bbox_inches='tight')