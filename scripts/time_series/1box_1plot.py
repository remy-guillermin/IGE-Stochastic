# LIBRAIRIES   
import numpy as np
import xarray as xr
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import cmocean
matplotlib.use('agg')

# PARAMETERS
# Path
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
simu = '/lus/store/CT1/c1601279/lweiss/run_croco/run_swio2_deter2_2017_202*'
file = 'swio_avg.nc'
obs = '/lus/scratch/CT1/c1601279/rguillermin/REGRIDDED/OBS'

# Boxes
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
cmap = cmocean.cm.tempo_r
n_colors = 4
colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
colors[0] = mcolors.to_rgba('palevioletred')
colors[-1] = mcolors.to_rgba('sandybrown')

fill_value = 9.96921e+36
roll = 9

SWIO = (25, 69, -36, 7)



fig, ax = plt.subplots(figsize=(6, 1))
ax.imshow([colors], extent=[0, n_colors, 0, 1])
ax.set_xticks([])
ax.set_yticks([])
plt.savefig(os.path.join(figures, f'colors.png'), dpi=300, transparent=True)
plt.close()



# FILES
croco = glob.glob(os.path.join(simu, file))
croco.sort()

sss_files = glob.glob(os.path.join(obs, 'SSS*'))
sss_files.sort()

sst_files = glob.glob(os.path.join(obs, 'SST*'))
sst_files.sort() 

zeta_files = glob.glob(os.path.join(obs, 'SLA*'))
zeta_files.sort()



# DATASET
obs_salt = xr.open_dataset(sss_files[0])["sos"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
    
obs_temp = xr.open_dataset(sst_files[0])["analysed_sst"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020'))) - 273.15
    
obs_zeta = xr.open_dataset(zeta_files[0])["sla"].sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))

obs_ds = xr.merge([obs_zeta, obs_temp, obs_salt]).rename({'x': 'xi_rho', 'y': 'eta_rho'}).isel(depth=0, drop=True)

croco_ds = xr.open_dataset(croco[0])[['temp', 'salt', 'zeta']].isel(s_rho=-1, drop=True)
croco_ds = croco_ds.where((croco_ds != fill_value), np.nan)



# GRID
g = xr.open_dataset(grid)[['lon_rho', 'lat_rho']]



for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
    print(name)
    region_mask = g.lon_rho.where((g.lon_rho > lon1) & (g.lon_rho < lon2) & (g.lat_rho > lat1) & (g.lat_rho < lat2), False)
    region_mask = region_mask.where((region_mask == False), True)
    
    region_obs = obs_ds.where(region_mask, drop=True).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    region_obs['time'] = region_obs.time.astype('datetime64[D]')
    region_obs = region_obs.groupby("time").reduce(np.nanmean)
    
    region_croco = croco_ds.where(region_mask, drop=True).mean(dim=['eta_rho', 'xi_rho'], skipna=True)
    
    # SALINITY
    fig, ax = plt.subplots(figsize=(10, 4))
    
    data = region_obs.sos
    time = region_obs.time
    
    ax.plot(time, data, 'r--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'r-', linewidth=2.5, label='OBS - MOG')
    
    data = region_croco.salt
    time = region_croco.time
    
    ax.plot(time, data, 'k--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'k-', linewidth=1.5, label='CROCO')
    
    fig.suptitle(f'SSS for {name} box')
    ax.set_ylabel('Salinity [psu]')
    ax.set_xlabel('Time')
    ax.tick_params("x", rotation=45)
    ax.set_xlim(region_croco.time.min(), region_croco.time.max())
    ax.grid(linewidth=0.3)
    ax.legend(loc='upper right')
    plt.savefig(os.path.join(figures, f'{name}_salinity.png'), dpi=300, transparent=True)
    plt.close()

    # TEMPERATURE
    fig, ax = plt.subplots(figsize=(10, 4))
    
    data = region_obs.analysed_sst
    time = region_obs.time
    
    ax.plot(time, data, 'r--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'r-', linewidth=2.5, label='OBS - OSTIA')
    
    data = region_croco.temp
    time = region_croco.time
    
    ax.plot(time, data, 'k--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'k-', linewidth=1.5, label='CROCO')
    
    fig.suptitle(f'SST for {name} box')
    ax.set_ylabel('Temperature [°C]')
    ax.set_xlabel('Time')
    ax.tick_params("x", rotation=45)
    ax.set_xlim(region_croco.time.min(), region_croco.time.max())
    ax.grid(linewidth=0.3)
    ax.legend(loc='lower right')
    plt.savefig(os.path.join(figures, f'{name}_temperature.png'), dpi=300, transparent=True)
    plt.close()

    # SALINITY ANOMALY
    fig, ax = plt.subplots(figsize=(10, 4))
    
    data = region_obs.sos - region_obs.sos.mean(dim='time')
    time = region_obs.time
    
    ax.plot(time, data, 'r--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'r-', linewidth=2.5, label='OBS - MOG')
    
    data = region_croco.salt - region_croco.salt.mean(dim='time')
    time = region_croco.time
    
    ax.plot(time, data, 'k--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'k-', linewidth=1.5, label='CROCO')
    
    fig.suptitle(f'SSA for {name} box')
    ax.set_ylabel('Surface Salinity Anomaly [psu]')
    ax.set_xlabel('Time')
    ax.tick_params("x", rotation=45)
    ax.set_xlim(region_croco.time.min(), region_croco.time.max())
    ax.grid(linewidth=0.3)
    ax.legend(loc='lower right')
    plt.savefig(os.path.join(figures, f'{name}_salinity_anomaly.png'), dpi=300, transparent=True)
    plt.close()

    # TEMPERATURE ANOMALY
    fig, ax = plt.subplots(figsize=(10, 4))
    
    data = region_obs.analysed_sst - region_obs.analysed_sst.mean(dim='time')
    time = region_obs.time
    
    ax.plot(time, data, 'r--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'r-', linewidth=2.5, label='OBS - OSTIA')
    
    data = region_croco.temp - region_croco.temp.mean(dim='time')
    time = region_croco.time
    
    ax.plot(time, data, 'k--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'k-', linewidth=1.5, label='CROCO')
    
    fig.suptitle(f'STA for {name} box')
    ax.set_ylabel('Surface Temperature Anomaly[°C]')
    ax.set_xlabel('Time')
    ax.tick_params("x", rotation=45)
    ax.set_xlim(region_croco.time.min(), region_croco.time.max())
    ax.grid(linewidth=0.3)
    ax.legend(loc='lower right')
    plt.savefig(os.path.join(figures, f'{name}_temperature_anomaly.png'), dpi=300, transparent=True)
    plt.close()

    # LEVEL ANOMALY
    fig, ax = plt.subplots(figsize=(10, 4))
    
    data = region_obs.analysed_sst - region_obs.analysed_sst.mean(dim='time')
    time = region_obs.time
    
    ax.plot(time, data, 'r--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'r-', linewidth=2.5, label='OBS - SSALTO')
    
    data = region_croco.temp - region_croco.temp.mean(dim='time')
    time = region_croco.time
    
    ax.plot(time, data, 'k--', linewidth=1)
    rolling_mean = np.convolve(data, np.ones(roll)/roll, mode='same')
    ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], 'k-', linewidth=1.5, label='CROCO')
    
    fig.suptitle(f'SLA for {name} box')
    ax.set_ylabel('Sea Level Anomaly [m]')
    ax.set_xlabel('Time')
    ax.tick_params("x", rotation=45)
    ax.set_xlim(region_croco.time.min(), region_croco.time.max())
    ax.grid(linewidth=0.3)
    ax.legend(loc='lower right')
    plt.savefig(os.path.join(figures, f'{name}_level_anomaly.png'), dpi=300, transparent=True)
    plt.close()
