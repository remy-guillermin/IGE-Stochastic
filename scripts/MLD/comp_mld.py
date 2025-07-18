# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import cartopy.crs as ccrs
import glob
import os
matplotlib.use('Agg')


# PARAMETERS
# Path
work = '/lus/work/CT1/c1601279/rguillermin'
sto_ens = ['run_swio2_stoens30_ini', 'run_swio2_stoens30_str']
grid = '/lus/work/CT1/c1601279/rguillermin/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'

lat_index = [430, 280]
SWIO = (25, 69, -36, 7)



# DATASET
ensemble_ini_files = sorted(glob.glob(os.path.join(work, 'MLD', sto_ens[0], '*')))
print(f"{len(ensemble_ini_files)} files found for {sto_ens[0]}.")
ensemble_str_files = sorted(glob.glob(os.path.join(work, 'MLD', sto_ens[1], '*')))
print(f"{len(ensemble_str_files)} files found for {sto_ens[1]}.")



ds_ini = xr.concat([xr.open_dataset(file) for file in ensemble_ini_files], dim='member')
print('ini concatenated.')
ds_str = xr.concat([xr.open_dataset(file) for file in ensemble_str_files], dim='member')
print('str concatenated.')



std_ini = ds_ini.std(dim='member')
std_str = ds_str.std(dim='member')
print('STD computed')

mean_ini = ds_ini.mean(dim='member')
mean_str = ds_str.mean(dim='member')
print('Mean computed')



g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

h = h.where((h != 50), 0)

# norm = mcolors.Normalize(vmin=0, vmax=0.3)

for i in range(min(ds_ini.time.size, ds_str.time.size)):
    pass
    # date = ds_ini.time.isel(time=i).data.astype('datetime64[D]')
    # print(f"Plotting {date}.")
    
    # fig, axs = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    # fig.suptitle(f"Ensemble INI standard deviation {date}")



date = ds_ini.time.isel(time=i).data.astype('datetime64[D]')
print(f"Plotting {date}.")

fig, axs = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"Ensemble standard deviation for the mixed layer depth - {date}")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands, zorder=3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}

# INI
ax = axs[0]
ax.set_title("INI")
norm = mcolors.Normalize(vmin=0, vmax=30)

pcm = ax.pcolormesh(lon, lat, std_ini.isel(time=i).mld, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree(), norm=norm)
cb = plt.colorbar(pcm, ax=ax, label='MLD standard deviation [m]', orientation='vertical')

# INI
ax = axs[1]
ax.set_title("STR")

pcm = ax.pcolormesh(lon, lat, std_str.isel(time=i).mld, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree(), norm=norm)
cb = plt.colorbar(pcm, ax=ax, label='MLD standard deviation [m]', orientation='vertical')


fig.tight_layout()
fig.savefig(os.path.join(figures, f'mld_std_comp_{date}'), dpi=300)
plt.close()



date = ds_ini.time.isel(time=i).data.astype('datetime64[D]')
print(f"Plotting {date}.")

fig, axs = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"Ensemble mean for the mixed layer depth - {date}")

for ax in axs:
    ax.set_extent(SWIO)
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

    land_color = ccrs.cartopy.feature.COLORS['land']

    minor_islands = ccrs.cartopy.feature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.add_feature(minor_islands, zorder=3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'color': 'k'}

# INI
ax = axs[0]
ax.set_title("INI")
norm = mcolors.Normalize(vmin=-70, vmax=0)

pcm = ax.pcolormesh(lon, lat, mean_ini.isel(time=i).mld, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree(), norm=norm)
cb = plt.colorbar(pcm, ax=ax, label='MLD [m]', orientation='vertical')

# INI
ax = axs[1]
ax.set_title("STR")

pcm = ax.pcolormesh(lon, lat, mean_str.isel(time=i).mld, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree(), norm=norm)
cb = plt.colorbar(pcm, ax=ax, label='MLD [m]', orientation='vertical')


fig.tight_layout()
fig.savefig(os.path.join(figures, f'mld_mean_comp_{date}'), dpi=300)
plt.close()
