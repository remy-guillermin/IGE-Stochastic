import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import cmocean
import sys
import matplotlib.colors as mcolors

# colors = ['sandybrown', 'goldenrod', 'seagreen', 'steelblue']
# colors = ["#4c72b0", "#dd8452", "#55a868", "#8172b3"]

cmap = cmocean.cm.tempo_r
n_colors = 4
colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
colors[0] = mcolors.to_rgba('palevioletred')
colors[-1] = mcolors.to_rgba('sandybrown')
fig, ax = plt.subplots(figsize=(6, 1))
ax.imshow([colors], extent=[0, n_colors, 0, 1])
ax.set_xticks([])
ax.set_yticks([])
plt.savefig(os.path.join('/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/run_swio2_stoens30_2017_ini', f'colors2.png'), dpi=100, bbox_inches='tight')
plt.close()

# sys.exit()

# Définition des chemins
store = '/lus/store/CT1/c1601279/lweiss/run_croco/'
scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE_dev/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'

stoens = 'run_swio2_stoens30_2017_ini'
# stoens = 'run_swio2_stoens30_20170101_01'
ref = 'run_swio2_deter2_2017_2023'
obs = 'OBS'

path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/' + stoens

# Chemins des fichiers d'ensemble et déterministe
ensemble_files = [os.path.join(scratch, stoens, f"{i:03d}swiose_avg.nc") for i in range(1, 31)]
deterministic_file = os.path.join(store, ref, 'swio_avg.nc')

# Noms des zones d'intérêt
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
# Ouverture des données
ensemble_datasets = [xr.open_dataset(f) for f in ensemble_files]
print('members : ', len(ensemble_datasets), 'time dim : ', len(ensemble_datasets[0]['time']))
deterministic_ds = xr.open_dataset(deterministic_file)
print('time deterministic dim : ', len(deterministic_ds['time']))

# Fonctions 
def region_subset(ds, var, box):
    "Calcule la moyenne spatiale de `var` dans `box` en fonction des dimensions spécifiques `lat_rho` et `lon_rho`."
    lon_min, lon_max, lat_min, lat_max = box
    if len(ds[var].dims) == 4:
        subset = ds[var].isel(s_rho=-1)  # Prend la couche de surface seulement si la variable est 3D
    else:
        subset = ds[var]

    subset = subset.where(
        (ds['lon_rho'] >= lon_min) & (ds['lon_rho'] <= lon_max) & 
        (ds['lat_rho'] >= lat_min) & (ds['lat_rho'] <= lat_max), drop=True)
    # spatial_mean = subset.mean(dim=['eta_rho', 'xi_rho'])
    # print(subset.shape) 
    return subset

def compute_variable_time_series(ds_list, var, box):
    "Calcule les séries temporelles pour chaque membre, la moyenne d'ensemble et la série déterministe pour `var` dans `box`."
    member_means = [region_subset(ds, var, box).mean(dim=['eta_rho', 'xi_rho']) for ds in ds_list]
    ## ensemble_mean2 = sum(member_means) / len(member_means) ## moyenne des moyennes spatiales 
    # Calculer la moyenne d'ensemble point à point
    ensemble_mean2d = sum([region_subset(ds, var, box) for ds in ds_list]) / len(ds_list)
    ensemble_mean = ensemble_mean2d.mean(dim=['eta_rho', 'xi_rho'])
    # deviation
    deviation2d = [region_subset(ds, var, box) - ensemble_mean2d for ds in ds_list]
    sq_deviation2d = [(dev ** 2) for dev in deviation2d]
    deviation = [dev.mean(dim=['eta_rho', 'xi_rho']) for dev in deviation2d]
    sq_deviation = [sq_dev.mean(dim=['eta_rho', 'xi_rho']) for sq_dev in sq_deviation2d]
    # RMS
    rms = np.sqrt(sum(sq_deviation) / len(ds_list))
    
    return member_means, ensemble_mean, deviation, rms 

# Variables et types de figures
variables = {'sst': 'temp', 'sss': 'salt', 'ssh': 'zeta'}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']
rms_list = {}
rms_uv_list = {}

# Boucle principale
for var_name, var in variables.items():
    print(var_name,var)
    for box, name in zip(boxes, names):
        print(box,name)
        # Calcul des séries temporelles
        member_means, ensemble_mean, deviation, rms = compute_variable_time_series(ensemble_datasets, var, box)
        
        # 1. Figure des moyennes
        plt.figure(figsize=(10, 4))
        for i, member in enumerate(member_means):
            plt.plot(member.time, member, label=f'Member {i+1}', color='grey', alpha=1, linewidth=0.3)
        plt.plot(ensemble_mean.time, ensemble_mean, color='indianred', label='Ensemble mean', linewidth=1)
        # plt.plot(deterministic_mean.time, deterministic_mean, 'r--', label='Deterministic', linewidth=1)
        plt.title(f'Average {var_name.upper()} in the {name} zone')
        # plt.legend()
        # plt.xlabel('Temps')
        plt.xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        plt.ylabel(var_name.upper())
        plt.grid(linewidth=0.3)
        plt.savefig(os.path.join(path_fig, f'{var_name}_{name}_mean_comparison.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # 2. Figure des écarts (x' = xi - x_mean)
        plt.figure(figsize=(10, 4))
        for i, member in enumerate(deviation):
            # deviation = member - ensemble_mean
            plt.plot(member.time, member, label=f'Member {i+1}', color='grey', alpha=1, linewidth=0.3)
        plt.title(f'{var_name.upper()} deviation from the ensemble mean in the {name} zone')
        # plt.xlabel('Temps')
        plt.xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        plt.ylabel(f'{var_name.upper()} deviation')
        # plt.legenid()
        plt.grid(linewidth=0.3)
        plt.savefig(os.path.join(path_fig, f'{var_name}_{name}_deviation.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # 3. Figure RMS des écarts        
        # rms = np.sqrt(sum([(member - ensemble_mean)**2 for member in member_means]) / len(member_means))
        rms_list[name] = []  # Initialise une liste vide pour cette clé
        rms_list[name].append(rms)

    # Tracer les courbes de RMS de toutes les zones sur la même figure    
    plt.figure(figsize=(10, 4))
    for i, (name, rms) in enumerate(rms_list.items()):
        plt.plot(rms[0].time, rms[0], label=f'{name}', color=colors[i], alpha=1)
    plt.title(f'RMS of {var_name.upper()} deviation from the ensemble mean')
    # plt.xlabel('Temps')
    plt.xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
    plt.ylabel(f'RMS of {var_name.upper()} deviation')
    plt.legend()
    plt.grid(linewidth=0.3)
    plt.savefig(os.path.join(path_fig, f'{var_name}_rms_deviation.png'), dpi=300, bbox_inches='tight')
    plt.close()

