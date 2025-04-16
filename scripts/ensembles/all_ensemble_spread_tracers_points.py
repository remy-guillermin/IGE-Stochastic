import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import cmocean
import glob
import matplotlib.colors as mcolors
import re
matplotlib.use('agg')



# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
work = '/lus/work/CT1/c1601279/rguillermin'
grid = os.path.join(work, 'grid/croco_grid_swio2.nc')
sto_ens = ['run_swio2_stoens30_201*_ini', 'run_swio2_stoens30_201*_CD']
obs = os.path.join(work, 'REGRIDDED', 'OBS')

ensembles = ['INI', 'STR']
linestyles = ['solid', 'dashed']

# Plot
points = [(290, 245), (285, 182)]
names = ['North-Mad', 'Islands']
cmap = cmocean.cm.tempo_r
n_colors = 2
box_colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
box_colors[0] = mcolors.to_rgba('palevioletred')
box_colors[-1] = mcolors.to_rgba('sandybrown')



# DATASET
ensemble_ini_files = sorted(glob.glob(os.path.join(work, 'NaN_CORRECTED', sto_ens[0], '*')))
print(f"{len(ensemble_ini_files)} files found for {sto_ens[0].replace('201*_', '')}.")
ensemble_str_files = sorted(glob.glob(os.path.join(work, 'NaN_CORRECTED', sto_ens[1], '*')))
print(f"{len(ensemble_str_files)} files found for {sto_ens[1].replace('201*_', '')}.")

pattern = re.compile(r'/(\d{3})swiose_avg.nc$')
ensemble_ini_dict = {}
ensemble_str_dict = {}

for path in ensemble_ini_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_ini_dict.keys():
            ensemble_ini_dict[key] = []
        ensemble_ini_dict[key].append(path)

for path in ensemble_str_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_str_dict.keys():
            ensemble_str_dict[key] = []
        ensemble_str_dict[key].append(path)



ensemble_ini_datasets = []       
for key, paths in ensemble_ini_dict.items():
    print(f"Merging member{key} for {sto_ens[0].replace('201*_', '')}.")
    ds = xr.Dataset()
    for path in paths:
        temp_ds = xr.open_dataset(path)
        mask = np.full(temp_ds.time.size, True)
        if ds.variables !=  {}:
            mask = ~np.isin(temp_ds.time, ds.time)
        ds = xr.merge([ds, temp_ds.isel(time=mask)])
    ensemble_ini_datasets.append(ds)

print('members : ', len(ensemble_ini_datasets), 'time dim : ', len(ensemble_ini_datasets[0]['time']))

ensemble_str_datasets = []       
for key, paths in ensemble_str_dict.items():
    print(f"Merging member{key} for {sto_ens[1].replace('201*_', '')}.")
    ds = xr.Dataset()
    for path in paths:
        temp_ds = xr.open_dataset(path)
        mask = np.full(temp_ds.time.size, True)
        if ds.variables !=  {}:
            mask = ~np.isin(temp_ds.time, ds.time)
        ds = xr.merge([ds, temp_ds.isel(time=mask)])
    ensemble_str_datasets.append(ds)

print('members : ', len(ensemble_str_datasets), 'time dim : ', len(ensemble_str_datasets[0]['time']))



sss_files = glob.glob(os.path.join(obs, 'SSS*'))
sss_files.sort()

sst_files = glob.glob(os.path.join(obs, 'SST*'))
sst_files.sort() 

zeta_files = glob.glob(os.path.join(obs, 'SLA*'))
zeta_files.sort()

# GRID
g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho']]

# DATASET
obs_salt = xr.open_dataset(sss_files[0]).sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
    
obs_temp = xr.open_dataset(sst_files[0]).sel(time = slice(np.datetime64('2017'), np.datetime64('2020'))) - 273.15
    
obs_zeta = xr.open_dataset(zeta_files[0]).sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
obs_ds = xr.merge([obs_zeta, obs_temp, obs_salt])

combined_ini = xr.concat(ensemble_ini_datasets, dim='ensemble')
combined_str = xr.concat(ensemble_str_datasets, dim='ensemble')



ensemble_mean2d_ini = combined_ini.mean(dim='ensemble')
ensemble_mean2d_str = combined_str.mean(dim='ensemble')
print('Ensemble mean computed.')



variables = {'sst': ('temp', 'analysed_sst'), 'sss': ('salt', 'sos'), 'ssh': ('zeta', 'adt')}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']



for (ilon, ilat), name in zip(points, names):        
    member_means_ini = combined_ini.isel(eta_rho=ilon, xi_rho=ilat)
    
    ensemble_mean_ini = ensemble_mean2d_ini.isel(eta_rho=ilon, xi_rho=ilat)
    
    member_deviations_ini = (member_means_ini - ensemble_mean_ini)
    sq_deviation_ini = member_deviations_ini ** 2
    
    member_means_str = combined_str.isel(eta_rho=ilon, xi_rho=ilat)

    ensemble_mean_str = ensemble_mean2d_str.isel(eta_rho=ilon, xi_rho=ilat)
    
    member_deviations_str = (member_means_str - ensemble_mean_str)
    sq_deviation_str = (member_deviations_str ** 2)
    
    obs_mean = obs_ds.isel(eta_rho=ilon, xi_rho=ilat)
    
    for var_name, (var, obs_var) in variables.items():
        print(name, var_name,var)
        
        # MEAN
        fig, ax = plt.subplots(figsize=(10, 5))
        
        for i in range(member_means_ini.sizes['ensemble']):
            member = member_means_ini.isel(ensemble=i)
            ax.plot(member.time, member[var], color='black', alpha=0.3, linewidth=0.3)
            
        ax.plot(member.time, member[var], color='black', alpha=0.3, label='INI - Member', linewidth=0.3)
        
        ax.plot(ensemble_mean_ini.time, ensemble_mean_ini[var], color='black', label='INI - Ensemble mean', linewidth=1)
        
        for i in range(member_means_str.sizes['ensemble']):
            member = member_means_str.isel(ensemble=i)
            ax.plot(member.time, member[var], color='royalblue', alpha=0.3, linewidth=0.3)
            
        ax.plot(member.time, member[var], color='royalblue', alpha=0.3, label='STR - Member', linewidth=0.3)
        
        ax.plot(ensemble_mean_str.time, ensemble_mean_str[var], color='royalblue', label='STR - Ensemble mean', linewidth=1)
        
        # ax.plot(obs_mean.time, obs_mean[obs_var], color='red', label='OBS', linewidth=2)
        
        fig.suptitle(f'Average {var_name.upper()} in the {name}')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_point_{var}.png'), dpi=300, transparent=True)
        plt.close()
        
        # DEVIATION
        fig, ax = plt.subplots(figsize=(10, 5))

        std_ini = member_deviations_str.std(dim='ensemble')
        ax.fill_between(std_ini.time, - std_ini[var], std_ini[var], color='royalblue', alpha=0.5, label='STR - Member', linewidth=0.5)

        std_ini = member_deviations_ini.std(dim='ensemble')
        ax.fill_between(std_ini.time, - std_ini[var], std_ini[var], color='firebrick', alpha=0.5, label='INI - Member', linewidth=0.5)
        
        fig.suptitle(f'{var_name.upper()} variance in the {name}')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
        ax.set_ylim(-np.max(np.abs(ax.get_ylim())), np.max(np.abs(ax.get_ylim())))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_point_{var}_deviation.png'), dpi=300, transparent=True)
        plt.close()
