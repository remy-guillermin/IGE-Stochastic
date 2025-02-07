import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cmocean
import cmcrameri

# Définir les chemins
scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swiose_0002'

# Définir les indices du point spécifique
lat_index = 249 # 249
lon_index = 195

# Charger les données
ds = xr.open_dataset(scratch + simu + '/swiose_his.nc')  # , engine='h5netcdf')
temp = ds['temp'][:, :, :, :]  # Température à tous les niveaux de profondeur au point spécifié
time = ds['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
s_rho = ds['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = ds['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = ds['hc'].values # S-coordinate parameter, critical depth
h = ds['h'][:, :]  # bathymetry at RHO-points
ds.close()

# Calcul z0 a nonlinear vertical transformation 
def calc_depth(s, Cs, hc, h):
    N = len(s_rho)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h ## (hc * s[k] + h * Cs[k])
    return depth

# Calcul de la profondeur
depth = calc_depth(s_rho, Cs_rho, hc, h)

# Assurer que les temps sont au format compréhensible si nécessaire (par exemple, en utilisant pandas)
# time = pd.to_datetime(time, unit='s')  # Ajustez le format de temps si nécessaire

# Convertir les données en numpy arrays pour faciliter le traitement
temp = temp[:, :, lat_index, lon_index].values
time = time.values
depth = depth[:, lat_index, lon_index]

cmap = cmocean.cm.thermal
# Créer la figure
fig, ax = plt.subplots(figsize=(10, 6))

# Créer une carte de chaleur
c = ax.pcolormesh(time, depth, temp.T, shading='auto', cmap=cmap)

# Ajouter une colorbar
cb = plt.colorbar(c, ax=ax)
cb.set_label('Temperature °C')

# Ajouter des labels et un titre
ax.set_xlabel('Time')
ax.set_ylabel('Depth')
ax.set_title(f'Vertical Temperature Profiles at lat index {lat_index}, lon index {lon_index}')

# Inverser l'axe des profondeurs si nécessaire (si les profondeurs augmentent avec l'index)
# ax.invert_yaxis()

# Afficher la figure
plt.grid(True)
plt.tight_layout()

# Sauvegarder la figure
plt.savefig(path_fig + '/temp_profiles_' + str(lat_index) + '_' + str(lon_index) + '_' + simu + '.png', dpi=300, bbox_inches='tight')

plt.show()