# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import os


# PARAMETERS
# Path
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
simu = '/lus/store/CT1/c1601279/lweiss/run_croco/run_swio2_deter2_2017_2023'
file = 'swio_avg.nc'

# Selection
date = '2018-06-15'
SWIO = (25, 69, -36, 7)

# Plot
figsize = (20, 6)
sst_cmap = cmocean.cm.thermal
sss_cmap = cmocean.cm.haline
ssh_cmap = cmcrameri.cm.roma_r

g = xr.open_dataset(grid)
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
g.close()
print("Grid loaded.")

d = xr.open_dataset(os.path.join(simu, file))
salt = d.salt
temp = d.temp
zeta = d.zeta
time = d.time
d.close()
print('Data loaded.')

if np.datetime64(date) not in time.values.astype('datetime64[D]'):
    print(f'{date} not found in simulation.')
    date = np.datetime_as_string(time.values.astype('datetime64[D]')[0])
    print(f'Setting {date} as the new date.')

fill_value = 9.96921e+36
salt = salt[:, -1, :, :]
salt = salt.where((salt != fill_value), np.nan)
zeta = zeta.where((zeta != fill_value), np.nan)
temp = temp[:, -1, :, :]
temp = temp.where((temp != fill_value), np.nan)
print('NaN values added')

salt_mean = salt.mean(dim='time')
zeta_mean = zeta.mean(dim='time')
temp_mean = temp.mean(dim='time')

salt = salt.sel(time=date).mean(dim='time')
zeta = zeta.sel(time=date).mean(dim='time')
temp = temp.sel(time=date).mean(dim='time')

print('Temporal mean computed.')

gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}

fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"SWIO {date}")

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

# Salinity
ax = axs[0]
ax.set_title("Salinity")

norm = mpl.colors.Normalize(vmin=34, vmax=36)
data = salt
cmap = sss_cmap
label = r'$\overline{\overline{SSS}}$ [psu]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
cb.set_ticks(np.arange(norm.vmin, norm.vmax + 0.5, 0.5))

# Free surface
ax = axs[1]
ax.set_title("Free Surface")

norm = mpl.colors.Normalize(vmin=0, vmax=1)
data = zeta
cmap = ssh_cmap
label = r'$\overline{\overline{SSH}}$ [m]'
pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Temperature
ax = axs[2]
ax.set_title("Temperature")

norm = mpl.colors.Normalize(vmin=20, vmax=30)
data = temp
cmap = sst_cmap
label = r'$\overline{\overline{SST}}$ [°C]'

pcm = ax.pcolormesh(lon[:, :], lat[:, :], data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')


fig.tight_layout()
filename = os.path.join(figures, f'physical_croco_{os.path.splitext(file)[0]}.png')
fig.savefig(filename, dpi=300, transparent=True)
print(f'{filename} saved.')
plt.show()
