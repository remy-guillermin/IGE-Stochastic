import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

# Définition des chemins
data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

stoens = 'run_swio_stoens03_2017'
ref = 'run_swio_ref2_2017'

path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/' + stoens

# Chemins des fichiers d'ensemble et déterministe
ensemble_files = [os.path.join(data, stoens, fname) for fname in ['001swiose_avg.nc', '002swiose_avg.nc', '003swiose_avg.nc']]
deterministic_file = os.path.join(data, ref, 'swio_avg.nc')

# Noms des zones d'intérêt
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'North-Moz', 'South-Moz', 'Mascarene']
colors = ['saddlebrown', 'goldenrod', 'seagreen', 'steelblue']
# Ouverture des données
ensemble_datasets = [xr.open_dataset(f) for f in ensemble_files]
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
    speed = np.sqrt(u_data**2 + v_data**2)  # Calcul de la norme de vitesse
    speed_da = xr.DataArray(speed, coords=ds['temp'][:, -1, :-1, :-1].coords, dims=ds['temp'][:, -1, :-1, :-1].dims)
    return region_subset(ds.assign(speed=speed_da), 'speed', box)

def compute_variable_time_series(ds_list, deterministic_ds, var, box):
    """Calcule les séries temporelles pour chaque membre, la moyenne d'ensemble et la série déterministe pour `var` dans `box`."""
    if var == 'speed':
        member_means = [compute_speed(ds, box).mean(dim=['eta_rho', 'xi_rho']) for ds in ds_list]
        ensemble_mean2d = sum([compute_speed(ds, box) for ds in ds_list]) / len(ds_list)
        ensemble_mean = ensemble_mean2d.mean(dim=['eta_rho', 'xi_rho'])
        deterministic_mean = compute_speed(deterministic_ds, box).mean(dim=['eta_rho', 'xi_rho'])
        deviation2d = [compute_speed(ds, box) - ensemble_mean2d for ds in ds_list]
        sq_deviation2d = [(dev ** 2) for dev in deviation2d]
        deviation = [dev.mean(dim=['eta_rho', 'xi_rho']) for dev in deviation2d]
        sq_deviation = [sq_dev.mean(dim=['eta_rho', 'xi_rho']) for sq_dev in sq_deviation2d]
        rms = np.sqrt(sum(sq_deviation) / len(ds_list))
    else:
        member_means = [region_subset(ds, var, box).mean(dim=['eta_rho', 'xi_rho']) for ds in ds_list]
        ## ensemble_mean2 = sum(member_means) / len(member_means) ## moyenne des moyennes spatiales 
        # Calculer la moyenne d'ensemble point à point
        ensemble_mean2d = sum([region_subset(ds, var, box) for ds in ds_list]) / len(ds_list)
        ensemble_mean = ensemble_mean2d.mean(dim=['eta_rho', 'xi_rho'])
        deterministic_mean = region_subset(deterministic_ds, var, box).mean(dim=['eta_rho', 'xi_rho'])
        # deviation
        deviation2d = [region_subset(ds, var, box) - ensemble_mean2d for ds in ds_list]
        sq_deviation2d = [(dev ** 2) for dev in deviation2d]
        deviation = [dev.mean(dim=['eta_rho', 'xi_rho']) for dev in deviation2d]
        sq_deviation = [sq_dev.mean(dim=['eta_rho', 'xi_rho']) for sq_dev in sq_deviation2d]
        # RMS
        rms = np.sqrt(sum(sq_deviation) / len(ds_list))
    return member_means, ensemble_mean, deterministic_mean, deviation, rms 

# Variables et types de figures
# variables = {'uv': 'speed'}
variables = {'uv': 'speed', 'sst': 'temp', 'sss': 'salt', 'ssh': 'zeta'} # , 'speed': 'speed'}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']
rms_list = {}
# Boucle principale
for var_name, var in variables.items():
    print(var_name,var)
    for box, name in zip(boxes, names):
        print(box,name)
        # Calcul des séries temporelles
        member_means, ensemble_mean, deterministic_mean, deviation, rms = compute_variable_time_series(ensemble_datasets, deterministic_ds, var, box)

        # 1. Figure des moyennes
        plt.figure(figsize=(8, 4))
        for i, member in enumerate(member_means):
            plt.plot(member.time, member, label=f'Member {i+1}', alpha=0.5)
        plt.plot(ensemble_mean.time, ensemble_mean, 'k-', label='Ensemble mean', linewidth=2)
        plt.plot(deterministic_mean.time, deterministic_mean, 'r--', label='Deterministic', linewidth=1)
        plt.title(f'Average {var_name.upper()} in the {name} zone')
        plt.legend()
        # plt.xlabel('Temps')
        plt.ylabel(var_name.upper())
        plt.grid()
        plt.savefig(os.path.join(path_fig, f'{var_name}_{name}_mean_comparison.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # 2. Figure des écarts (x' = xi - x_mean)
        plt.figure(figsize=(8, 4))
        for i, member in enumerate(deviation):
            # deviation = member - ensemble_mean
            plt.plot(member.time, member, label=f'Member {i+1}', alpha=0.5)
        plt.title(f'{var_name.upper()} deviation from the ensemble mean in the {name} zone')
        # plt.xlabel('Temps')
        plt.ylabel(f'{var_name.upper()} deviation')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(path_fig, f'{var_name}_{name}_deviation.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # 3. Figure RMS des écarts        
        # rms = np.sqrt(sum([(member - ensemble_mean)**2 for member in member_means]) / len(member_means))
        rms_list[name] = []  # Initialise une liste vide pour cette clé
        rms_list[name].append(rms)

    # Tracer les courbes de RMS de toutes les zones sur la même figure    
    plt.figure(figsize=(8, 4))
    for i, (name, rms) in enumerate(rms_list.items()):
        plt.plot(rms[0].time, rms[0], label=f'{name}', color=colors[i], alpha=0.7)
    plt.title(f'RMS of {var_name.upper()} deviation from the ensemble mean')
    # plt.xlabel('Temps')
    plt.ylabel(f'RMS of {var_name.upper()} deviation')
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(path_fig, f'{var_name}_rms_deviation.png'), dpi=300, bbox_inches='tight')
    plt.close()