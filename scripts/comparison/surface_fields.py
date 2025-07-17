# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
import cmocean
import cartopy.crs as ccrs
import glob
import os
from functools import partial

# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures/Ensembles'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
obs_dir = '/lus/store/CT1/c1601279/rguillermin/REGRIDDED/OBS'
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
data = '/lus/store/CT1/c1601279/lweiss/run_croco/SWIO/run_swio2_deter2_2017_2023'

ensembles = ['INI', 'STR']
linestyles = ['solid', 'dashed']

date = '2017-06-15'
trimonth = {"DJF": slice('2016-12', '2017-02'), 'MAM': slice('2017-03', '2017-05'), 'JJA': slice('2017-06', '2017-08'), 'SON': slice('2017-09', '2017-11')}
months = 'JJA'
year = '2019'

SWIO = (25, 69, -36, 7)
# Plot
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}
figsize = (15, 6)
cmap = cmocean.cm.balance



data_files = glob.glob(os.path.join(data, f'swio_avg_{year}.nc'))
data_files.sort()

sss_files = glob.glob(os.path.join(obs_dir, 'SSS*'))
sss_files.sort()

sst_files = glob.glob(os.path.join(obs_dir, 'SST*'))
sst_files.sort()



obs_salt = xr.open_dataset(sss_files[0]).sel(time = slice(np.datetime64('2017'), np.datetime64('2020')))
    
obs_temp = xr.open_dataset(sst_files[0]).sel(time = slice(np.datetime64('2017'), np.datetime64('2020'))) - 273.15
    
obs_ds = xr.merge([obs_temp, obs_salt])

print("Observation datasets loaded.")

data_ds = xr.open_dataset(data_files[0])[['temp', 'salt']].isel(s_rho=-1)


data_ds['time'] = data_ds.time.astype('datetime64[D]')
obs_ds['time'] = obs_ds.time.astype('datetime64[D]')


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

data_ds = data_ds.where(msk)

fig, axs = plt.subplots(1, 2, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"Deviation from Observation {year}")

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

member = data_ds.sel(time=year)
obs = obs_ds.sel(time=year)[['analysed_sst', "sos"]].rename({'analysed_sst': 'temp', 'sos': 'salt'})

diff = (member - obs).mean(dim='time')
maximum = np.abs(diff).max()

# Salinity
ax = axs[0]
ax.set_title("Salinity")

pcm = ax.pcolormesh(lon, lat, diff.salt, cmap=cmap, transform=ccrs.PlateCarree(), vmin=-maximum.salt, vmax=maximum.salt)
cb = plt.colorbar(pcm, ax=ax, label='Salinity deviation [psu]', orientation='vertical')

# Salinity
ax = axs[1]
ax.set_title("Temperature")

pcm = ax.pcolormesh(lon, lat, diff.temp, cmap=cmap, transform=ccrs.PlateCarree(), vmin=-maximum.temp, vmax=maximum.temp)
cb = plt.colorbar(pcm, ax=ax, label='Temperature deviation [°C]', orientation='vertical')

fig.tight_layout()
filename = os.path.join(figures, f"deviation_obs_{year}_deter.png")
fig.savefig(filename, dpi=300)
print(f'{filename} saved.')
plt.show()
