import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import cmocean
import glob
import matplotlib
import matplotlib.colors as mcolors
import re
matplotlib.use('agg')



# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
work = '/lus/work/CT1/c1601279/rguillermin/'
nan_rem = '/lus/work/CT1/c1601279/rguillermin/NaN_CORRECTED/'
stoens = ['run_swio2_stoens30_201*_ini', 'run_swio2_stoens30_201*_str', 'run_swio2_stoens30_201*_gls']
obs = '/lus/store/CT1/c1601279/rguillermin/REGRIDDED/OBS'
grid = os.path.join(work, 'grid/croco_grid_swio2.nc')

ensembles = ['INI', 'STR', 'GLS']
linestyles = ['solid', 'dashed', 'dotted']

# Plot
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
cmap = cmocean.cm.tempo_r
n_colors = 4
box_colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
box_colors[0] = mcolors.to_rgba('palevioletred')
box_colors[-1] = mcolors.to_rgba('sandybrown')



ensemble_ini_files = sorted(glob.glob(os.path.join(nan_rem, stoens[0], "*swiose_avg.nc")))
print(f"{len(ensemble_ini_files)} files found for {stoens[0].replace('201*_', '')}.")
ensemble_str_files = sorted(glob.glob(os.path.join(nan_rem, stoens[1], "*swiose_avg.nc")))
print(f"{len(ensemble_str_files)} files found for {stoens[1].replace('201*_', '')}.")
ensemble_gls_files = sorted(glob.glob(os.path.join(nan_rem, stoens[2], "*swiose_avg.nc")))
print(f"{len(ensemble_gls_files)} files found for {stoens[2].replace('201*_', '')}.")

# GRID
g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho']]

pattern = re.compile(r'/(\d{3})swiose_avg.nc$')
ensemble_ini_dict = {}
ensemble_str_dict = {}
ensemble_gls_dict = {}

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

for path in ensemble_gls_files:
    match = pattern.search(path)
    if match:
        key = match.group(1)
        if key not in ensemble_gls_dict.keys():
            ensemble_gls_dict[key] = []
        ensemble_gls_dict[key].append(path)


ensemble_ini_datasets = []       
for key, paths in ensemble_ini_dict.items():
    print(f"Merging member{key} for {stoens[0].replace('201*_', '')}.")
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
    print(f"Merging member{key} for {stoens[1].replace('201*_', '')}.")
    ds = xr.Dataset()
    for path in paths:
        temp_ds = xr.open_dataset(path)
        mask = np.full(temp_ds.time.size, True)
        if ds.variables !=  {}:
            mask = ~np.isin(temp_ds.time, ds.time)
        ds = xr.merge([ds, temp_ds.isel(time=mask)])
    ensemble_str_datasets.append(ds)

print('members : ', len(ensemble_str_datasets), 'time dim : ', len(ensemble_str_datasets[0]['time']))

ensemble_gls_datasets = []       
for key, paths in ensemble_gls_dict.items():
    print(f"Merging member{key} for {stoens[2].replace('201*_', '')}.")
    ds = xr.Dataset()
    for path in paths:
        temp_ds = xr.open_dataset(path)
        mask = np.full(temp_ds.time.size, True)
        if ds.variables !=  {}:
            mask = ~np.isin(temp_ds.time, ds.time)
        ds = xr.merge([ds, temp_ds.isel(time=mask)])
    ensemble_gls_datasets.append(ds)

print('members : ', len(ensemble_gls_datasets), 'time dim : ', len(ensemble_gls_datasets[0]['time']))

combined_ini = xr.concat(ensemble_ini_datasets, dim='ensemble')
combined_str = xr.concat(ensemble_str_datasets, dim='ensemble')
combined_gls = xr.concat(ensemble_gls_datasets, dim='ensemble')

ensemble_mean2d_ini = combined_ini.mean(dim='ensemble')
ensemble_mean2d_str = combined_str.mean(dim='ensemble')
ensemble_mean2d_gls = combined_gls.mean(dim='ensemble')
print('Ensemble mean computed.')

variables = {'sst': 'temp', 'sss': 'salt'}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']
rms_list = {ensemble : {} for ensemble in ensembles}

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
    
    rms_ini = np.sqrt(sq_deviation_ini.sum(dim='ensemble') / len(sq_deviation_ini))
    rms_list['INI'][name] = rms_ini
    
    member_means2d_str = combined_str.where(region_mask, drop=True)
    member_means_str = member_means2d_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)

    ensemble_mean_region_str = ensemble_mean2d_str.where(region_mask, drop=True)
    ensemble_mean_str = ensemble_mean_region_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    member_deviations2d_str = (member_means2d_str - ensemble_mean_region_str)
    member_deviations_str = member_deviations2d_str.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    sq_deviation_str = (member_deviations2d_str ** 2).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    rms_str = np.sqrt(sq_deviation_str.sum(dim='ensemble') / len(sq_deviation_str))
    rms_list['STR'][name] = rms_str
    
    member_means2d_gls = combined_gls.where(region_mask, drop=True)
    member_means_gls = member_means2d_gls.mean(dim=['eta_rho', 'xi_rho'], skipna=True)

    ensemble_mean_region_gls = ensemble_mean2d_gls.where(region_mask, drop=True)
    ensemble_mean_gls = ensemble_mean_region_gls.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    member_deviations2d_gls = (member_means2d_gls - ensemble_mean_region_gls)
    member_deviations_gls = member_deviations2d_gls.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    sq_deviation_gls = (member_deviations2d_gls ** 2).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    rms_gls = np.sqrt(sq_deviation_gls.sum(dim='ensemble') / len(sq_deviation_gls))
    rms_list['GLS'][name] = rms_gls
       
    
    for var_name, var, in variables.items():
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
        
        for i in range(member_means_gls.sizes['ensemble']):
            member = member_means_gls.isel(ensemble=i)
            ax.plot(member.time, member[var], color='crimson', alpha=0.3, linewidth=0.3)
            
        ax.plot(member.time, member[var], color='crimson', alpha=0.3, label='GLS - Member', linewidth=0.3)
        
        ax.plot(ensemble_mean_gls.time, ensemble_mean_gls[var], color='crimson', label='GLS - Ensemble mean', linewidth=1)
        
        
        fig.suptitle(f'Average {var_name.upper()} in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_{var}.png'), dpi=300, transparent=True)
        plt.close()
        
        # DEVIATION
        fig, ax = plt.subplots(figsize=(10, 5))

        std_gls = member_deviations_gls.std(dim='ensemble')
        ax.fill_between(std_gls.time, - std_gls[var], std_gls[var], color='crimson', alpha=0.5, label='GLS', linewidth=1)
        
        std_str = member_deviations_str.std(dim='ensemble')
        ax.fill_between(std_str.time, - std_str[var], std_str[var], color='royalblue', alpha=0.5, label='STR', linewidth=1)

        std_ini = member_deviations_ini.std(dim='ensemble')
        ax.fill_between(std_ini.time, - std_ini[var], std_ini[var], color='black', alpha=0.5, label='INI', linewidth=1)
        
        fig.suptitle(f'{var_name.upper()} deviation from the ensemble mean in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
        if var == 'salt':
            ax.set_ylim(-0.06, 0.06)
        elif var == 'temp':
            ax.set_ylim(-0.2, 0.2)
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='upper right')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'all_ens_{name}_{var}_deviation.png'), dpi=300, transparent=True)
        plt.close()

for var_name, var in variables.items():
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (ensemble, ensemble_rms_list) in enumerate(rms_list.items()):
        for j, (name, rms) in enumerate(ensemble_rms_list.items()):
            ax.plot(rms[var].time, rms[var], label=name if (ensemble == 'INI') else None, linestyle=linestyles[i], color=box_colors[j])
            
            
    fig.suptitle(f'RMS of {var_name.upper()} deviation from the ensemble mean')
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2019-12-31'))
    ax.set_ylabel(f'RMS of {var_name.upper()} deviation')
    ax.tick_params("x", rotation=45)
    ax.legend(loc='upper left')
    ax.grid(linewidth=0.3)
    fig.tight_layout()
    plt.savefig(os.path.join(figures, f'all_ens_{var_name}_rms.png'), dpi=300, transparent=True)
    plt.close()






