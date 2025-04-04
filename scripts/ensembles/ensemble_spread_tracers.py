import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import cmocean
import glob
import sys
import matplotlib.colors as mcolors


# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
store = '/lus/store/CT1/c1601279/lweiss/run_croco/'
scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE_dev/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
stoens = 'run_swio2_stoens30_2017_ini'
ref = 'run_swio2_deter2_2017_2023'
obs = '/lus/store/CT1/c1601279/rguillermin/REGRIDDED/OBS'
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

# Plot
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
cmap = cmocean.cm.tempo_r
n_colors = 4
colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
colors[0] = mcolors.to_rgba('palevioletred')
colors[-1] = mcolors.to_rgba('sandybrown')

fill_value = 9.96921e+36
roll = 9



ensemble_files = [os.path.join(scratch, stoens, f"{i:03d}swiose_avg.nc") for i in range(1, 31)]
deterministic_file = os.path.join(store, ref, 'swio_avg_2017.nc')

sss_files = glob.glob(os.path.join(obs, 'SSS*'))
sss_files.sort()

sst_files = glob.glob(os.path.join(obs, 'SST*'))
sst_files.sort() 

zeta_files = glob.glob(os.path.join(obs, 'SLA*'))
zeta_files.sort()

# GRID
g = xr.open_dataset(grid)[['lon_rho', 'lat_rho']]



ensemble_datasets = []
for f in ensemble_files:
    print(f'Converting NaN in {f}')
    ds = xr.open_dataset(f)[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
    ds = ds.where((ds != fill_value), np.nan)
    ensemble_datasets.append(ds)

print('members : ', len(ensemble_datasets), 'time dim : ', len(ensemble_datasets[0]['time']))
deterministic_ds = xr.open_dataset(deterministic_file)[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
deterministic_ds = deterministic_ds.where((deterministic_ds != fill_value), np.nan)
print('time deterministic dim : ', len(deterministic_ds['time']))



# DATASET
obs_salt = xr.open_dataset(sss_files[0])["sos"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
    
obs_temp = xr.open_dataset(sst_files[0])["analysed_sst"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020'))) - 273.15
    
obs_zeta = xr.open_dataset(zeta_files[0])["sla"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
obs_ds = xr.merge([obs_zeta, obs_temp, obs_salt])



variables = {'sst': ('temp', 'analysed_sst'), 'sss': ('salt', 'sos'), 'ssh': ('zeta', 'adt')}
fig_types = ['mean_comparison', 'deviation', 'rms_deviation']
rms_list = {}
rms_uv_list = {}

for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
    print((lon1, lon2, lat1, lat2), name)
        
    region_mask = g.lon_rho.where((g.lon_rho > lon1) & (g.lon_rho < lon2) & (g.lat_rho > lat1) & (g.lat_rho < lat2), False)
    region_mask = region_mask.where((region_mask == False), True)
    
    member_means2d = [ds.where(region_mask, drop=True) for ds in ensemble_datasets]
    member_means = [member.mean(dim=['eta_rho', 'xi_rho'], skipna=True) for member in member_means2d]
    ensemble_mean2d = sum(member_means2d) / len(member_means)
    ensemble_mean = ensemble_mean2d.mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    member_deviations2d = [ds.where(region_mask, drop=True) - ensemble_mean2d for ds in ensemble_datasets]
    member_deviations = [member.mean(dim=['eta_rho', 'xi_rho'], skipna=True) for member in member_deviations2d]
    
    sq_deviation2d = [(dev ** 2) for dev in member_deviations2d]
    sq_deviation = [member.mean(dim=['eta_rho', 'xi_rho'], skipna=True) for member in sq_deviation2d]
    
    rms = np.sqrt(sum(sq_deviation) / len(sq_deviation))
    rms_list[name] = rms
    
    obs_mean = obs_ds.where(region_mask, drop=True).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    for var_name, (var, obs_var) in variables.items():
        print(var_name,var)
        
        # MEAN
        fig, ax = plt.subplots(figsize=(10, 4))
        
        for member in member_means:
            ax.plot(member.time, member[var], color='grey', linewidth=0.3)
            
        ax.plot(member.time, member[var], color='grey', label='Member', linewidth=0.3)
        
        ax.plot(ensemble_mean.time, ensemble_mean[var], color='black', label='Ensemble mean', linewidth=1)
        
        if var_name != 'ssh':
            ax.plot(obs_mean.time, obs_mean[obs_var], color='red', label='OBS', linewidth=2)
        
        
        
        fig.suptitle(f'Average {var_name.upper()} in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='best')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'{stoens}_{name}_{var}.png'), dpi=300, transparent=True)
        plt.close()
        
        # DEVIATION     
        fig, ax = plt.subplots(figsize=(10, 4))

        for member in member_deviations:
            ax.plot(member.time, member[var], color='grey', linewidth=0.3)

        ax.plot(member.time, member[var], color='grey', label='Member', alpha=1, linewidth=0.3)
        
        
        fig.suptitle(f'{var_name.upper()} deviation from the ensemble mean in the {name} zone')        
        ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
        ax.set_ylim(-np.max(np.abs(ax.get_ylim())), np.max(np.abs(ax.get_ylim())))
        ax.set_ylabel(var_name.upper())
        ax.tick_params("x", rotation=45)
        ax.grid(linewidth=0.3)
        ax.legend(loc='best')
        fig.tight_layout()
        plt.savefig(os.path.join(figures, f'{stoens}_{name}_{var}_deviation.png'), dpi=300, transparent=True)
        plt.close()


for var_name, (var, obs_var) in variables.items():
    fig, ax = plt.subplots(figsize=(10, 4))

    for i, (name, rms) in enumerate(rms_list.items()):
        ax.plot(rms[var].time, rms[var], label=f'{name}', color=colors[i])

    fig.suptitle(f'RMS of {var_name.upper()} deviation from the ensemble mean')
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2017-12-31'))
    ax.set_ylabel(f'RMS of {var_name.upper()} deviation')
    ax.tick_params("x", rotation=45)
    ax.legend(loc='upper left')
    ax.grid(linewidth=0.3)
    fig.tight_layout()
    plt.savefig(os.path.join(figures, f'{stoens}_{var_name}_rms.png'), dpi=300, transparent=True)
    plt.close()
