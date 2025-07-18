# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import glob
import os
from matplotlib.gridspec import GridSpec
matplotlib.use('Agg')

work = '/lus/work/CT1/c1601279/rguillermin'
store = '/lus/store/CT1/c1601279/rguillermin'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'

stoens = 'run_croco/run_swio2_stoens30_gls_2017'

grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'

lat_index = [430, 280]
SWIO = (25, 69, -36, 7)
CRT1 = 0.03 # kg/m3


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

ensemble_path = glob.glob(os.path.join(store, stoens, '*avg.nc'))
print(f"{len(ensemble_path)} files found for {stoens}.")
ds = xr.open_dataset(ensemble_path[0])

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

ds = xr.open_dataset(ensemble_path[0])
depth_sigma = calc_depth(ds.s_rho, ds.Cs_rho, ds.hc, h)
s_rho = ds.s_rho
temp = ds.temp
ds.close()

depth_level = xr.Dataset(
    data_vars=dict(
        depth_level=(['s_rho', 'eta_rho', 'xi_rho'], depth_sigma),
    ),
    coords=dict(
        eta_rho=eta_rho,
        xi_rho=xi_rho,
        s_rho=s_rho,
    )
)

h = h.where((h != 50), 0)

rho = ds.rho.where((g.mask_rho), np.nan)

MLDC = rho - rho.isel(s_rho=-2) < CRT1 # Condition on density for MLD

fig = plt.figure(figsize=(10,6))
gs = GridSpec(2, 2, width_ratios=[20,1], height_ratios=[1,1], wspace=0.05)

date = MLDC.isel(time=-1).time.data.astype('datetime64[D]')

rho_ref = r'\rho_{ref}'

fig.suptitle(rf'Zonal slices for the density anomaly criterion $\rho - {rho_ref} < {CRT1}$ $kg/m³$ - {date}')

ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
axs = [ax0, ax1]

land_color = ccrs.cartopy.feature.COLORS['land']

for i, (ax, id) in enumerate(zip(axs, lat_index)):
    
    Lon = lon.isel(eta_rho=id, xi_rho=slice(0,431))
    Lon = Lon.expand_dims({'s_rho':50})
    Lon['s_rho'] = ds.s_rho
    
    levels = depth_level.depth_level.isel(eta_rho=id, xi_rho=slice(0,431))
    
    data = MLDC.sel(eta_rho=id+1).isel(time=-1, xi_rho=slice(0,431))
    pcm = ax.pcolormesh(Lon, levels, data)
        
    ax.set_xlim((40,65))
    ax.fill_between(lon[id, :], -h[id, :].values, y2=min(-h[id, :]), color=land_color)
    ax.plot(lon[id, :], -h[id, :].values, color='black', linewidth=0.5)
    
    ax.set_ylim(-100, 0)
    
    ticks = ax.get_yticks()[1:-1]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
    ax.set_title(f'Slice S{i+1}')

plt.setp(ax0.get_xticklabels(), visible=False)
    
ticks = ax1.get_xticks()
ax1.set_xticks(ticks)
ax1.set_xticklabels([f'{tick:.0f}°E' for tick in ticks])
ax1.set_xlabel('Longitude')

fig.tight_layout()
fig.savefig(f'{figures}/slices_mask_mld_{date}.png', dpi=300, transparent=True)
plt.close()


mask = - MLDC.sum(dim='s_rho')
mask['eta_rho'] = mask.eta_rho - 1
mask['xi_rho'] = mask.xi_rho - 1
# mask = mask.where((mask !=0), np.nan)
MLD = depth_level.isel(s_rho=mask)
MLD = MLD.where((g.mask_rho), np.nan)



fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

fig.suptitle(f"Mixing density layer - {date}")

ax.set_extent(SWIO)   
ax.coastlines(resolution='50m')
ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'k'}
gl.ylabel_style = {'size': 10, 'color': 'k'}

data = MLD.isel(time=-1).depth_level

pcm = ax.pcolormesh(lon, lat, data, cmap=cmocean.cm.deep, transform=ccrs.PlateCarree())
cb = fig.colorbar(pcm, ax=ax, label='Mixing layer depth [m]')

fig.tight_layout()
fig.savefig(f'{figures}/mld_{date}.png')
plt.close()


fig, axs = plt.subplots(2, 1, figsize=(10,6))

fig.suptitle(f'Zonal slices for the mixing layer - {date}')

land_color = ccrs.cartopy.feature.COLORS['land']

for i, (ax, id) in enumerate(zip(axs, lat_index)):
    
    data = MLD.depth_level.sel(eta_rho=id).isel(time=-1)
    
    ax.plot(lon[id, :], data)
        
    ax.set_xlim((40,65))
    ax.fill_between(lon[id, :], -h[id, :].values, y2=min(-h[id, :]), color=land_color)
    ax.plot(lon[id, :], -h[id, :].values, color='black', linewidth=0.5)
    
    ax.set_ylim(-75, 0)
    
    ticks = ax.get_yticks()[1:-1]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick)} m' for tick in ticks])
    
    ax.set_title(f'Slice S{i+1}')

plt.setp(ax0.get_xticklabels(), visible=False)
    
ticks = ax1.get_xticks()
ax1.set_xticks(ticks)
ax1.set_xticklabels([f'{tick:.0f}°E' for tick in ticks])
ax1.set_xlabel('Longitude')

fig.tight_layout()
fig.savefig(f'{figures}/slices_mld_{date}.png', dpi=300, transparent=True)
plt.close()
