# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import glob
import os
import re
matplotlib.use('agg')



# PATH
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles/Frames'
work = '/lus/work/CT1/c1601279/rguillermin'
grid = os.path.join(work, 'grid/croco_grid_swio2.nc')
sto_ens = ['run_swio2_stoens30_201*_ini', 'run_swio2_stoens30_201*_CD']

# PLOT OPTIONS
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}
figsize = (15, 6)
cmap = cmocean.cm.dense
SWIO = (25, 69, -36, 7)

# Create Directories
os.makedirs(figures, exist_ok=True)

for ens in sto_ens:
    path = ens.replace('201*_', '')
    os.makedirs(os.path.join(figures, path), exist_ok=True)
    


# Search and open files
g = xr.open_dataset(grid)
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
pm = g['pm'][:,:] 
pn = g['pn'][:,:]
msk_inv = np.where(msk == 0, msk, np.nan)
h = g['h'][:, :]
angle = g['angle'][:, :]
g.close()
print("Grid loaded.")



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
for _, paths in ensemble_ini_dict.items():
    ds = xr.open_mfdataset(paths)
    ensemble_ini_datasets.append(ds)

print('members : ', len(ensemble_ini_datasets), 'time dim : ', len(ensemble_ini_datasets[0]['time']))

ensemble_str_datasets = []       
for _, paths in ensemble_str_dict.items():
    ds = xr.Dataset()
    for path in paths:
        ds = xr.merge([ds, xr.open_dataset(path)], compat='override')
    ensemble_str_datasets.append(ds)

print('members : ', len(ensemble_str_datasets), 'time dim : ', len(ensemble_str_datasets[0]['time']))



combined_ini = xr.concat(ensemble_ini_datasets, dim='ensemble')
combined_str = xr.concat(ensemble_str_datasets, dim='ensemble')



for i in range(combined_ini.time.size):
    date = combined_ini.time.isel(time=i).data.astype('datetime64[D]')
    print(f"Plotting {date}.")
    
    fig, axs = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(f"Ensemble INI standard deviation {date}")

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
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'color': 'k'}

    ensemble = combined_ini.isel(time=i)
    std = ensemble.std(dim='ensemble')
    
    # Salinity
    ax = axs[0]
    ax.set_title("Salinity")

    pcm = ax.pcolormesh(lon, lat, std.salt, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree())
    cb = plt.colorbar(pcm, ax=ax, label='Salinity standard deviation [psu]', orientation='vertical')

    # Free surface
    ax = axs[1]
    ax.set_title("Free Surface")

    pcm = ax.pcolormesh(lon, lat, std.zeta, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree())
    cb = plt.colorbar(pcm, ax=ax, label='Free surface standard deviation [m]', orientation='vertical')

    # Temperature
    ax = axs[2]
    ax.set_title("Temperature")

    pcm = ax.pcolormesh(lon, lat, std.temp, cmap=cmocean.cm.dense, transform=ccrs.PlateCarree())
    cb = plt.colorbar(pcm, ax=ax, label='Temperature standard deviation [°C]', orientation='vertical')

    fig.tight_layout()
    filename = os.path.join(figures, sto_ens[0].replace('201*_', ''), f"deviation_{date}_ensemble.png")
    fig.savefig(filename, dpi=300, transparent=True)
    plt.close()