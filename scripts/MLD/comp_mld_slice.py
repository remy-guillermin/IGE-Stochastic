
# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import os
matplotlib.use('Agg')



# PARAMETERS
# Path
work = '/lus/work/CT1/c1601279/rguillermin'
sto_ens = ['run_swio2_stoens30_ini', 'run_swio2_stoens30_str', 'run_swio2_stoens30_gls']
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
depth_lvl = '/lus/work/CT1/c1601279/rguillermin/grid/croco_depth_level.nc'

lat_index = [430, 280]
SWIO = (25, 69, -36, 7)

colors = ['black', 'royalblue', 'crimson']



# DATASET
ensemble_ini_files = sorted(glob.glob(os.path.join(work, 'MLD', sto_ens[0], '*')))
print(f"{len(ensemble_ini_files)} files found for {sto_ens[0]}.")
ensemble_str_files = sorted(glob.glob(os.path.join(work, 'MLD', sto_ens[1], '*')))
print(f"{len(ensemble_str_files)} files found for {sto_ens[1]}.")
ensemble_gls_files = sorted(glob.glob(os.path.join(work, 'MLD', sto_ens[2], '*')))
print(f"{len(ensemble_gls_files)} files found for {sto_ens[2]}.")



def calc_depth(
    s: xr.DataArray, 
    Cs: xr.DataArray, 
    hc: float, 
    h: xr.DataArray
) -> xr.DataArray:
    """
    Compute depth using the S-coordinate transformation.

    Parameters
    ----------
    s : xr.DataArray
        S-coordinate at RHO-points, typically ranging from -1 (surface) to 0 (bottom).
    Cs : xr.DataArray
        S-coordinate stretching curves at W-points, defining the vertical stretching.
    hc : float
        Critical depth parameter (in meters), influencing vertical terrain-following transformation.
    h : xr.DataArray
        Bathymetric depth at RHO-points (in meters), representing the seafloor depth.

    Returns
    -------
    xr.DataArray
        Computed depth at RHO-points.
    """
    N = len(s)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h
    return depth



g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

h = h.where((h != 50), 0)

depth_level = xr.open_dataset(depth_lvl)
print(f'depth level file opened: {depth_lvl}')



ds_ini = [xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_ini_files]
ds_str = [xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_str_files]
ds_gls = [xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_gls_files]



index = 250 
# index = min(ds_ini[0].time.size, ds_str[0].time.size, ds_gls[0].time.size) - 1

fig, axs = plt.subplots(3, 2, figsize=(16,9), sharex=True, sharey=True)

date = ds_ini[0].isel(time=index).time.data.astype('datetime64[D]')

fig.suptitle(rf'Zonal slice MLD - {date}')

land_color = ccrs.cartopy.feature.COLORS['land']

axs[0,0].set_title('Slice S1')
axs[0,1].set_title('Slice S2')

for ax in axs[:,0]:        
    ax.fill_between(lon[lat_index[0], :], -h[lat_index[0], :].values, y2=min(-h[lat_index[0], :]), color=land_color)
    ax.plot(lon[lat_index[0], :], -h[lat_index[0], :].values, color='black', linewidth=0.5)
    
    ax.set_xlim((40,65))
    ax.set_ylim(-100, 0)
    
for ax in axs[:,1]:        
    ax.fill_between(lon[lat_index[1], :], -h[lat_index[1], :].values, y2=min(-h[lat_index[1], :]), color=land_color)
    ax.plot(lon[lat_index[1], :], -h[lat_index[1], :].values, color='black', linewidth=0.5)
    
    ax.set_xlim((40,65))
    
    
ini_mean_mld = []
for i in range(len(lat_index)):
    mean_mld = xr.zeros_like(ds_ini[0].mld.isel(time=index, eta_rho=i))
    for ds in ds_ini:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld += mld
    mean_mld /= len(ds_ini)
    ini_mean_mld.append(mean_mld)
        
ini_std_mld = [[], []]
for i, ax in enumerate(axs[0,:]):
    id = lat_index[i]
    for ds in ds_ini:
        mld = ds.mld.isel(time=index, eta_rho=i)
        std_mld = mld - ini_mean_mld[i]
        ini_std_mld[i].append(std_mld.data)
        
        ax.plot(lon[id, :], mld, color=colors[0], linewidth=0.3, alpha=0.3)
        
    ini_std_mld[i] = np.concatenate(ini_std_mld[i])
    rms = np.sqrt(np.nansum((ini_std_mld[i] ** 2)) / (~np.isnan(ini_std_mld[i])).sum())
    ax.plot(lon[id, :], ini_mean_mld[i], color=colors[0], label=f'INI - RMSE: {rms:.1f} m')
  
  
str_mean_mld = []
for i in range(len(lat_index)):
    mean_mld = xr.zeros_like(ds_str[0].mld.isel(time=index, eta_rho=i))
    for ds in ds_str:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld += mld
    mean_mld /= len(ds_str)
    str_mean_mld.append(mean_mld)
 
str_std_mld = [[], []]   
for i, ax in enumerate(axs[1,:]):
    id = lat_index[i]
    for ds in ds_str:
        mld = ds.mld.isel(time=index, eta_rho=i)
        std_mld = mld - str_mean_mld[i]
        str_std_mld[i].append(std_mld.data)
        ax.plot(lon[id, :], mld, color=colors[1], linewidth=0.3, alpha=0.3)
        
    str_std_mld[i] = np.concatenate(str_std_mld[i])
    rms = np.sqrt(np.nansum((str_std_mld[i] ** 2)) / (~np.isnan(str_std_mld[i])).sum())
    ax.plot(lon[id, :], str_mean_mld[i], color=colors[1], label=f'STR - RMSE: {rms:.1f} m')
    
      
gls_mean_mld = []
for i in range(len(lat_index)):
    mean_mld = xr.zeros_like(ds_gls[0].mld.isel(time=index, eta_rho=i))
    for ds in ds_gls:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld += mld
    mean_mld /= len(ds_gls)
    gls_mean_mld.append(mean_mld)

gls_std_mld = [[], []]
for i, ax in enumerate(axs[2,:]):
    id = lat_index[i]
    for ds in ds_gls:
        mld = ds.mld.isel(time=index, eta_rho=i)
        std_mld = mld - gls_mean_mld[i]
        gls_std_mld[i].append(std_mld.data)
        ax.plot(lon[id, :], mld, color=colors[2], linewidth=0.3, alpha=0.3)
        
    gls_std_mld[i] = np.concatenate(gls_std_mld[i])
    rms = np.sqrt(np.nansum((gls_std_mld[i] ** 2)) / (~np.isnan(gls_std_mld[i])).sum())
    ax.plot(lon[id, :], gls_mean_mld[i], color=colors[2], label=f'GLS - RMSE: {rms:.1f} m')

for ax in axs[:,0]:
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
for ax in axs.flatten():
    ax.grid(alpha=0.5)
    ax.legend(loc='best')
    
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])

fig.tight_layout()
fig.savefig(f'{figures}/slices_mld_comp_{date}.png', dpi=300)
plt.close()


fig, axs = plt.subplots(1, 2, figsize=(16,4), sharex=True, sharey=True)

date = ds_ini[0].isel(time=index).time.data.astype('datetime64[D]')

fig.suptitle(f'Zonal slice MLD - {date}')

land_color = ccrs.cartopy.feature.COLORS['land']

ax = axs[0]
ax.set_title('Slice S1')
ax.fill_between(lon[lat_index[0], :], -h[lat_index[0], :].values, y2=min(-h[lat_index[0], :]), color=land_color)
ax.plot(lon[lat_index[0], :], -h[lat_index[0], :].values, color='black', linewidth=0.5)

ax = axs[1]
ax.set_title('Slice S1')      
ax.fill_between(lon[lat_index[1], :], -h[lat_index[1], :].values, y2=min(-h[lat_index[1], :]), color=land_color)
ax.plot(lon[lat_index[1], :], -h[lat_index[1], :].values, color='black', linewidth=0.5)

ini_mean_mld = []
str_mean_mld = []
gls_mean_mld = []
for i in range(2):
    id = lat_index[i]
    mean_mld_ini = xr.zeros_like(ds_ini[0].mld.isel(time=index, eta_rho=i))
    mean_mld_str = xr.zeros_like(ds_str[0].mld.isel(time=index, eta_rho=i))
    mean_mld_gls = xr.zeros_like(ds_gls[0].mld.isel(time=index, eta_rho=i))
    
    for ds in ds_ini:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld_ini += mld
        
    mean_mld_ini /= 30
    ini_mean_mld.append(mean_mld_ini)
    
    for ds in ds_str:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld_str += mld
        
    mean_mld_str /= 30
    str_mean_mld.append(mean_mld_str)
    
    for ds in ds_gls:
        mld = ds.mld.isel(time=index, eta_rho=i)
        mean_mld_gls += mld
        
    mean_mld_gls /= 30
    gls_mean_mld.append(mean_mld_gls)


for i, (ax, mld) in enumerate(zip(axs, ini_mean_mld)):
    id = lat_index[i]
    ax.plot(lon[id, :], mld, color=colors[0], label='INI')

for i, (ax, mld) in enumerate(zip(axs, str_mean_mld)):
    id = lat_index[i]
    ax.plot(lon[id, :], mld, color=colors[1], label='STR')

for i, (ax, mld) in enumerate(zip(axs, gls_mean_mld)):
    id = lat_index[i]
    ax.plot(lon[id, :], mld, color=colors[2], label='GLS')

for ax in axs:
    ax.set_xlim((40,65))
    ax.set_ylim(-100, 0)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
    ticks = ax.get_xticks()
    
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])
    
    ax.grid(alpha=0.5)
    ax.legend(loc='lower right')



fig.tight_layout()
fig.savefig(f'{figures}/slices_mld_comp_sup_{date}.png', dpi=300)
plt.close()


fig, axs = plt.subplots(1, 2, figsize=(16,4))

fig.suptitle('Bathymetry')

land_color = ccrs.cartopy.feature.COLORS['land']

ax = axs[0]
ax.set_title('Slice S1')
ax.fill_between(lon[lat_index[0], :], -h[lat_index[0], :].values, y2=-6000, color=land_color)
ax.plot(lon[lat_index[0], :], -h[lat_index[0], :].values, color='black', linewidth=0.5)


ax = axs[1]
ax.set_title('Slice S1')      
ax.fill_between(lon[lat_index[1], :], -h[lat_index[1], :].values, y2=-6000, color=land_color)
ax.plot(lon[lat_index[1], :], -h[lat_index[1], :].values, color='black', linewidth=0.5)

for ax in axs:
    ax.set_xlim((40,65))
    ax.set_ylim((-5500, 0))
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
    ticks = ax.get_xticks()
    
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(tick)}°E' for tick in ticks])
    
    ax.grid(alpha=0.5)
    
fig.tight_layout()
fig.savefig(f'{figures}/slices_bathy.png', dpi=300)
plt.close()