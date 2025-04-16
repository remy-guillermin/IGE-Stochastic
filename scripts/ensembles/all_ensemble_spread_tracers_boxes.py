import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import cmocean
import glob
import sys
import matplotlib.colors as mcolors
import time



# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
nan_rem = '/lus/work/CT1/c1601279/rguillermin/NaN_CORRECTED/'
stoens = ['run_swio2_stoens30_2017_ini', 'run_swio2_stoens30_2017_CD']
obs = '/lus/store/CT1/c1601279/rguillermin/REGRIDDED/OBS'
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

ensembles = ['INI', 'STR']
linestyles = ['solid', 'dashed']

# Plot
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
cmap = cmocean.cm.tempo_r
n_colors = 4
box_colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
box_colors[0] = mcolors.to_rgba('palevioletred')
box_colors[-1] = mcolors.to_rgba('sandybrown')



ensemble_ini_files = [os.path.join(nan_rem, stoens[0], f"{i:03d}swiose_avg.nc") for i in range(1, 31)]
ensemble_str_files = [os.path.join(nan_rem, stoens[1], f"{i:03d}swiose_avg.nc") for i in range(1, 31)]

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



ensemble_ini_datasets = []
start = time.time()
for f in ensemble_ini_files:
    start = time.time()
    ds = xr.open_dataset(f)
    ensemble_ini_datasets.append(ds)

print('members : ', len(ensemble_ini_datasets), 'time dim : ', len(ensemble_ini_datasets[0]['time']))

ensemble_str_datasets = []
for f in ensemble_str_files:
    start = time.time()
    ds = xr.open_dataset(f)
    ensemble_str_datasets.append(ds)

print('members : ', len(ensemble_str_datasets), 'time dim : ', len(ensemble_str_datasets[0]['time']))



combined_ini = xr.concat(ensemble_ini_datasets, dim='ensemble')
combined_str = xr.concat(ensemble_str_datasets, dim='ensemble')



ensemble_mean2d_ini = combined_ini.mean(dim='ensemble')
ensemble_mean2d_str = combined_str.mean(dim='ensemble')
print('Ensemble mean computed.')

# member_deviations2d_ini = combined_ini - ensemble_mean2d_ini
# member_deviations2d_str = combined_ini - ensemble_mean2d_str
# print('Member deviation computed.')

# sq_deviation2d_ini = member_deviations2d_ini ** 2
# sq_deviation2d_str = member_deviations2d_str ** 2
# print('Member STD computed.')



variables = {'sst': ('temp', 'analysed_sst'), 'sss': ('salt', 'sos'), 'ssh': ('zeta', 'adt')}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']
rms_list = {ensemble : {} for ensemble in ensembles}
rms_uv_list = {}



for (lon1, lon2, lat1, lat2), name in zip(boxes, names):        
    region_mask = g.lon_rho.where((g.lon_rho > lon1) & (g.lon_rho < lon2) & (g.lat_rho > lat1) & (g.lat_rho < lat2), False)
    region_mask = region_mask.where((region_mask == False), True)
    
    member_means2d_ini = combined_ini.where(region_mask, drop=True)
    member_means_ini = member_means2d_ini.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    ensemble_mean_region_ini = ensemble_mean2d_ini.where(region_mask, drop=True)
    ensemble_mean_ini = ensemble_mean_region_ini.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    member_deviations2d_ini = (member_means2d_ini - ensemble_mean_region_ini)
    member_deviations_ini = member_deviations2d_ini.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    sq_deviation_ini = (member_deviations2d_ini ** 2).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    rms_ini = sq_deviation_ini.sum(dim='ensemble') / len(sq_deviation_ini)
    rms_list['INI'][name] = rms_ini
    
    member_means2d_str = combined_str.where(region_mask, drop=True)
    member_means_str = member_means2d_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)

    ensemble_mean_region_str = ensemble_mean2d_str.where(region_mask, drop=True)
    ensemble_mean_str = ensemble_mean_region_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    member_deviations2d_str = (member_means2d_str - ensemble_mean_region_str)
    member_deviations_str = member_deviations2d_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    sq_deviation_str = (member_deviations2d_str ** 2).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    rms_str = sq_deviation_str.sum(dim='ensemble') / len(sq_deviation_str)
    rms_list['STR'][name] = rms_str
       
    obs_mean = obs_ds.where(region_mask, drop=True).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    for var_name, (var, obs_var) in variables.items():
        print(name, var_name,var)
        
        # MEAN
        fig, ax = plt.subplots(figsize=(10, 4))
        
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
        
        ax.plot(obs_mean.time, obs_mean[obs_var], color='red', label='OBS', linewidth=2)
        
        fig.suptitle(f'Average {var_name.upper()} in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_{var}.png'), dpi=300, transparent=True)
        plt.close()
        
        # DEVIATION
        fig, ax = plt.subplots(figsize=(10, 4))

        std_ini = member_deviations_ini.std(dim='ensemble')
        ax.fill_between(std_ini.time, - std_ini[var], std_ini[var], color='firebrick', alpha=0.5, label='INI - Member', linewidth=0.5)

        std_ini = member_deviations_str.std(dim='ensemble')
        ax.fill_between(std_ini.time, - std_ini[var], std_ini[var], color='royalblue', alpha=0.5, label='STR - Member', linewidth=0.5)
        
        fig.suptitle(f'{var_name.upper()} deviation from the ensemble mean in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        ax.set_ylim(-np.max(np.abs(ax.get_ylim())), np.max(np.abs(ax.get_ylim())))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_{var}_deviation.png'), dpi=300, transparent=True)
        plt.close()

for var_name, (var, obs_var) in variables.items():
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (ensemble, ensemble_rms_list) in enumerate(rms_list.items()):
        for j, (name, rms) in enumerate(ensemble_rms_list.items()):
            ax.plot(rms[var].time, rms[var], label=name if (ensemble == 'INI') else None, linestyle=linestyles[i], color=box_colors[j])
            
            
    fig.suptitle(f'RMS of {var_name.upper()} deviation from the ensemble mean')
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
    ax.set_ylabel(f'RMS of {var_name.upper()} deviation')
    ax.tick_params("x", rotation=45)
    ax.legend(loc='upper left')
    ax.grid(linewidth=0.3)
    fig.tight_layout()
    plt.savefig(os.path.join(figures, f'all_ens_{var_name}_rms.png'), dpi=300, transparent=True)
    plt.close()






