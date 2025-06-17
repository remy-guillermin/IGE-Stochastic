# LIBRAIRIES
import matplotlib.gridspec
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmocean
import cartopy.crs as ccrs
import glob
import os
from functools import partial
%matplotlib inline



# PARAMETERS
# Path
work = '/lus/work/CT1/c1601279/rguillermin'
sto_ens = ['run_swio2_stoens30_ini', 'run_swio2_stoens30_CD', 'run_swio2_stoens30_gls']
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'
sim = '/lus/work/CT1/c1601279/rguillermin/grid/swio_avg.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'

lat_index = [430, 280]
points = [(280, 380), (280, 190), (430, 170)]
points_names = ['Mascarene', 'Islands', 'Somali']
SWIO = (25, 69, -36, 7)

colors = ['black', 'royalblue', 'crimson']
transparent = False



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
Lon = g.lon_rho
Lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

ds = xr.open_dataset(sim)[['s_rho', 'Cs_rho', 'hc', 'eta_rho', 'xi_rho']]
depth_sigma = calc_depth(ds.s_rho, ds.Cs_rho, ds.hc, h)

depth_level = xr.Dataset(
    data_vars=dict(
        depth_level=(['s_rho', 'eta_rho', 'xi_rho'], depth_sigma),
    ),
    coords=dict(
        eta_rho=ds.eta_rho,
        xi_rho=ds.xi_rho,
        s_rho=ds.s_rho,
    )
)

h = h.where((h != 50), 0)
ds.close()



ds_ini = xr.concat([xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_ini_files], dim='member')
print('INI concatenated')
ds_str = xr.concat([xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_str_files], dim='member')
print('STR concatenated')
ds_gls = xr.concat([xr.open_dataset(file).isel(eta_rho=lat_index) for file in ensemble_gls_files], dim='member')
print('GLS concatenated')



mean_ini = ds_ini.mean(dim='member')
std_ini = ds_ini.std(dim='member')
print('INI computed')

mean_str = ds_str.mean(dim='member')
std_str = ds_str.std(dim='member')
print('STR computed')

mean_gls = ds_gls.mean(dim='member')
std_gls = ds_gls.std(dim='member')
print('GLS computed')



slices_mean_ini = mean_ini.mean(dim='xi_rho')
slices_mean_str = mean_str.mean(dim='xi_rho')
slices_mean_gls = mean_gls.mean(dim='xi_rho')

slices_var_ini = (std_ini ** 2).mean(dim='xi_rho')
slices_var_str = (std_str ** 2).mean(dim='xi_rho')
slices_var_gls = (std_gls ** 2).mean(dim='xi_rho')

slices_rmse_ini = np.sqrt(slices_var_ini.mean(dim='time'))
slices_rmse_str = np.sqrt(slices_var_str.mean(dim='time'))
slices_rmse_gls = np.sqrt(slices_var_gls.mean(dim='time'))

slices_mean_ini_ensemble = ds_ini.mean(dim='xi_rho')
slices_mean_str_ensemble = ds_str.mean(dim='xi_rho')
slices_mean_gls_ensemble = ds_gls.mean(dim='xi_rho')

max_date = min(ds_ini.time.max(), ds_str.time.max(), ds_gls.time.max())

for i, id in enumerate(lat_index):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f'Mean MLD for Slice S{i+1}')
    
    ini = slices_mean_ini.isel(eta_rho=i)
    members_ini = slices_mean_ini_ensemble.isel(eta_rho=i)
    
    str = slices_mean_str.isel(eta_rho=i)
    members_str = slices_mean_str_ensemble.isel(eta_rho=i)
    
    gls = slices_mean_gls.isel(eta_rho=i)
    members_gls = slices_mean_gls_ensemble.isel(eta_rho=i, time=slice(0,330))
    
    
    for j in range(members_ini.member.size):
        member = members_ini.isel(member=j)
        ax.plot(member.time, member.mld, color=colors[0], linewidth=0.3, alpha=0.3)
        
    ax.plot(ini.time, ini.mld, color=colors[0], label=f'INI - RMSE: {slices_rmse_ini.mld.isel(eta_rho=i):.1f} m')
    
    
    for j in range(members_str.member.size):
        member = members_str.isel(member=j)
        ax.plot(member.time, member.mld, color=colors[1], linewidth=0.3, alpha=0.3)
        
    ax.plot(str.time, str.mld, color=colors[1], label=f'STR - RMSE: {slices_rmse_str.mld.isel(eta_rho=i):.1f} m')
    
    
    for j in range(members_gls.member.size):
        member = members_gls.isel(member=j)
        ax.plot(member.time, member.mld, color=colors[2], linewidth=0.3, alpha=0.3)
        
    ax.plot(gls.time, gls.mld, color=colors[2], label=f'GLS - RMSE: {slices_rmse_gls.mld.isel(eta_rho=i):.1f} m')
        
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2018-01-01'))
    ax.set_xlim(np.datetime64('2017-01-01'), max_date)
    ax.set_ylim(-80, 0)
    
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    ax.tick_params("x", rotation=45)
    
    ax.legend(loc='best')
    ax.grid(alpha=0.5)
    fig.tight_layout()
    fig.savefig(f'{figures}/slice_S{i+1}_mld_comp_sup_time_series.png', dpi=300, transparent=transparent)
    plt.show()
    
    
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f'Deviation of MLD for Slice S{i+1}')
    
    ini_deviation = np.sqrt(slices_var_ini.isel(eta_rho=i))
    str_deviation = np.sqrt(slices_var_str.isel(eta_rho=i))
    gls_deviation = np.sqrt(slices_var_gls.isel(eta_rho=i))
    
    ax.fill_between(gls_deviation.time, - gls_deviation.mld, gls_deviation.mld, color=colors[2], alpha=0.5, label=f'GLS - RMSE: {slices_rmse_gls.mld.isel(eta_rho=i):.1f} m', linewidth=1)
    ax.fill_between(str_deviation.time, - str_deviation.mld, str_deviation.mld, color=colors[1], alpha=0.5, label=f'STR - RMSE: {slices_rmse_str.mld.isel(eta_rho=i):.1f} m', linewidth=1)
    ax.fill_between(ini_deviation.time, - ini_deviation.mld, ini_deviation.mld, color=colors[0], alpha=0.5, label=f'INI - RMSE: {slices_rmse_ini.mld.isel(eta_rho=i):.1f} m', linewidth=1)
    
    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2018-01-01'))
    ax.set_xlim(np.datetime64('2017-01-01'), max_date)
    ax.set_ylim(-20, 20)
    
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{int(tick)} m' for tick in ticks])
    ax.tick_params("x", rotation=45)
    
    ax.legend(loc='best')
    ax.grid(alpha=0.5)
    fig.tight_layout()
    fig.savefig(f'{figures}/slice_S{i+1}_mld_deviation_comp_sup_time_series.png', dpi=300, transparent=transparent)
    plt.show()



max_date = min(ds_ini.time.max(), ds_str.time.max(), ds_gls.time.max())

for i, ((lon_id, lat_id), name) in enumerate(zip(points, points_names)):
    print(i, lon_id, lat_id, name)
    lon, lat = Lon[lon_id, lat_id].data, Lat[lon_id, lat_id].data
    print(f'lon: {lon} - lat: {lat}')
    
    ini = ds_ini.sel(eta_rho=lon_id+1, xi_rho=lat_id+1)
    str = ds_str.sel(eta_rho=lon_id+1, xi_rho=lat_id+1)
    gls = ds_gls.sel(eta_rho=lon_id+1, xi_rho=lat_id+1)
    
    ini_mean = ini.mean(dim='member')
    ini_std = ini.std(dim='member')
    
    str_mean = str.mean(dim='member')
    str_std = str.std(dim='member')
    
    gls_mean = gls.mean(dim='member')
    gls_std = gls.std(dim='member')
    
    fig = plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
    axs = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1], projection=ccrs.PlateCarree())]
    fig.suptitle(f'Mean MLD for {name} Point')
    
    ax = axs[0]
        
    ax.plot(ini_mean.time, ini_mean.mld, color=colors[0], label=f"INI - RMSE: {np.sqrt((ini_std.mld**2).mean(dim='time').data):.1f} m")
        
    ax.plot(str_mean.time, str_mean.mld, color=colors[1], label=f"STR - RMSE: {np.sqrt((str_std.mld**2).mean(dim='time').data):.1f} m")

    ax.plot(gls_mean.time, gls_mean.mld, color=colors[2], label=f"GLS - RMSE: {np.sqrt((gls_std.mld**2).mean(dim='time').data):.1f} m")

    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2018-01-01'))
    ax.set_xlim(np.datetime64('2017-01-01'), max_date)
    ax.set_ylim(-100, 0)
    
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    ax.tick_params("x", rotation=45)
    
    ax.legend(loc='best')
    ax.grid(alpha=0.5)
    
    ax = axs[1]
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black')
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.OCEAN)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands)
    
    ax.scatter(lon, lat, c='red')
    
    ax.set_extent(SWIO)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='-', linewidth=0.2, color='k', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'k'}
    gl.ylabel_style = {'size': 10, 'color': 'k'}

    fig.tight_layout()
    fig.savefig(f'{figures}/points_{name}_mld_comp_sup_time_series.png', dpi=300, transparent=transparent)
    plt.show()

    fig = plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
    axs = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1], projection=ccrs.PlateCarree())]
    fig.suptitle(f'Mean MLD for {name} Point')
    
    ax = axs[0]
    
    ax.fill_between(gls_std.time, - gls_std.mld, gls_std.mld, color=colors[2], alpha=0.5, label=f"GLS - RMSE: {np.sqrt((gls_std.mld**2).mean(dim='time').data):.1f} m", linewidth=1)
    
    ax.fill_between(str_std.time, - str_std.mld, str_std.mld, color=colors[1], alpha=0.5, label=f"STR - RMSE: {np.sqrt((str_std.mld**2).mean(dim='time').data):.1f} m", linewidth=1)
    
    ax.fill_between(ini_std.time, - ini_std.mld, ini_std.mld, color=colors[0], alpha=0.5, label=f"INI - RMSE: {np.sqrt((ini_std.mld**2).mean(dim='time').data):.1f} m", linewidth=1)

    ax.set_xlim(np.datetime64('2017-01-01'), np.datetime64('2018-01-01'))
    ax.set_xlim(np.datetime64('2017-01-01'), max_date)
    ax.set_ylim(-30, 30)
    
    ticks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    ax.tick_params("x", rotation=45)
    
    ax.legend(loc='best')
    ax.grid(alpha=0.5)
    
    ax = axs[1]
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black')
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.OCEAN)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands)
    
    ax.scatter(lon, lat, c='red')
    
    ax.set_extent(SWIO)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='-', linewidth=0.2, color='k', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'k'}
    gl.ylabel_style = {'size': 10, 'color': 'k'}

    fig.tight_layout()
    fig.savefig(f'{figures}/points_{name}_mld_deviation_comp_sup_time_series.png', dpi=300, transparent=transparent)
    plt.show()