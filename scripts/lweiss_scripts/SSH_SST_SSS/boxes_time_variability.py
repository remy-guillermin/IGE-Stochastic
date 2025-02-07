import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

# Définition des chemins
data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

ref = 'run_swio_ref_2017_2023'

path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/' + ref

# Chemins des fichiers déterministe
deterministic_file = os.path.join(data, ref, 'swio_avg.nc')

# Noms des zones d'intérêt
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)]
names = ['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
# colors = ['saddlebrown', 'goldenrod', 'seagreen', 'steelblue']
colors = ['sandybrown', 'plum', 'cornflowerblue', 'mediumaquamarine']
# Ouverture des données
deterministic_ds = xr.open_dataset(deterministic_file)

# Fonctions utilitaires
def region_subset(ds, var, box):
    """Calcule la moyenne spatiale de `var` dans `box` en fonction des dimensions spécifiques `lat_rho` et `lon_rho`."""
    lon_min, lon_max, lat_min, lat_max = box
    if len(ds[var].dims) == 4:
        subset = ds[var].isel(s_rho=-1)  # Prend la couche de surface seulement si la variable est 3D
    else:
        subset = ds[var]

    subset = subset.where(
        (ds['lon_rho'] >= lon_min) & (ds['lon_rho'] <= lon_max) & 
        (ds['lat_rho'] >= lat_min) & (ds['lat_rho'] <= lat_max), drop=True
    )
    # spatial_mean = subset.mean(dim=['eta_rho', 'xi_rho']) 
    return subset

def compute_speed(ds, box):
    """Calcule la vitesse résultante sqrt(u^2 + v^2) en ajustant les dimensions de grille pour u et v."""
    u_data = ds['u'][:, -1, :-1, :].data  # Ajustement pour `u`
    v_data = ds['v'][:, -1, :, :-1].data  # Ajustement pour `v`
    speed = np.sqrt(u_data ** 2 + v_data ** 2)  # Calcul de la norme de vitesse
    speed_da = xr.DataArray(speed, coords=ds['temp'][:, -1, :-1, :-1].coords, dims=ds['temp'][:, -1, :-1, :-1].dims)
    return region_subset(ds.assign(speed=speed_da), 'speed', box)

def compute_ke(ds, box):
    u_data = ds['u'][:, -1, :-1, :].data  # Ajustement pour `u`
    v_data = ds['v'][:, -1, :, :-1].data  # Ajustement pour `v`
    ke = 0.5 * (u_data ** 2 + v_data ** 2) ## * area2[:,:-1, :-1].values
    ke_da = xr.DataArray(ke, coords=ds['temp'][:, -1, :-1, :-1].coords, dims=ds['temp'][:, -1, :-1, :-1].dims)
    return region_subset(ds.assign(ke=ke_da), 'ke', box)

def compute_variable_time_series(deterministic_ds, var, box):
    """Calcule les séries temporelles pour chaque membre, la moyenne d'ensemble et la série déterministe pour `var` dans `box`."""
    if var == 'speed':
        deterministic_mean = compute_speed(deterministic_ds, box).mean(dim=['eta_rho', 'xi_rho'])
    elif var == 'ke':
        deterministic_mean = compute_ke(deterministic_ds, box).mean(dim=['eta_rho', 'xi_rho'])
    elif var == 'sla':
        zeta_subset = region_subset(deterministic_ds, 'zeta', box)
        zeta_time_mean = zeta_subset.mean(dim='time')
        sla = zeta_subset - zeta_time_mean
        deterministic_mean = sla.mean(dim=['eta_rho', 'xi_rho'])
    else:
        # Calculer la moyenne d'ensemble point à point
        deterministic_mean = region_subset(deterministic_ds, var, box).mean(dim=['eta_rho', 'xi_rho'])
    return deterministic_mean 

# Variables et types de figures
# variables = {'uv': 'speed'}
variables = {'sla': 'sla', 'uv': 'speed', 'ke': 'ke', 'sst': 'temp', 'sss': 'salt', 'ssh': 'zeta'} # , 'speed': 'speed'}

data_list = {}
# Boucle principale
for var_name, var in variables.items():
    print(var_name,var)
    for box, name in zip(boxes, names):
        print(box,name)
        # Calcul des séries temporelles
        deterministic_mean = compute_variable_time_series(deterministic_ds, var, box)

        # 1. Figure des moyennes
        plt.figure(figsize=(8, 4))
        plt.plot(deterministic_mean.time, deterministic_mean, 'k', label='', linewidth=1)
        plt.fill_between(deterministic_mean.time, deterministic_mean-deterministic_mean.std(), deterministic_mean+deterministic_mean.std(), color='k', alpha=0.1)
        plt.title(f'Average {var_name.upper()} in the {name} zone')
        plt.legend(loc=0, ncol=1, fontsize=9)
        plt.xticks(fontsize=10)
        # plt.xlabel('Temps')
        plt.ylabel(var_name.upper(), fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.8)
        plt.savefig(os.path.join(path_fig, f'{var_name}_{name}_time_mean.png'), dpi=300, bbox_inches='tight')
        plt.close()

        data_list[name] = []
        data_list[name].append(deterministic_mean)        

    # Tracer les courbes de toutes les zones sur la même figure    
    plt.figure(figsize=(8, 4))
    for i, (name, data) in enumerate(data_list.items()):
        plt.plot(data[0].time, data[0], label=f'{name}', color=colors[i], alpha=0.8, linewidth=1.2)
        plt.fill_between(data[0].time, data[0]-data[0].std(), data[0]+data[0].std(), color=colors[i], alpha=0.2)
    # plt.xlabel('Temps')
    plt.xticks(fontsize=10)
    plt.ylabel(f'{var_name.upper()}', fontsize=14)
    plt.legend(loc=0, ncol=1, fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.8)
    plt.savefig(os.path.join(path_fig, f'{var_name}_time_mean.png'), dpi=300, bbox_inches='tight')
    plt.close()


