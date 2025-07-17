# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import glob
import os
from functools import partial



# PARAMETERS
# Path
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
data = '/lus/work/CT1/c1601279/rguillermin/RAW'

# Selection
date = '2018-06-15'
SWIO = (25, 69, -36, 7)

# Plot
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.5}
figsize = (20, 6)
sst_cmap = cmocean.cm.thermal
sss_cmap = cmocean.cm.haline
ssh_cmap = cmcrameri.cm.roma_r

year = np.datetime64(date).astype('datetime64[Y]')

sss_files = glob.glob(data + f'/SSS*/*{year}*')
sst_files = glob.glob(data + f'/SST*/*{year}*')
sla_files = glob.glob(data + f'/SLA*/*{year}*')

print(f"SSS files: {len(sss_files)}")
print(f"SST files: {len(sst_files)}")
print(f"SLA files: {len(sla_files)}")

def _preprocess(x, lon_bnds, lat_bnds):
    return x.sel(lon=slice(*lon_bnds), lat=slice(*lat_bnds))

lon_bnds, lat_bnds = (SWIO[0]-5, SWIO[1]+5), (SWIO[2]-5, SWIO[3]+5)
partial_func = partial(_preprocess, lon_bnds=lon_bnds, lat_bnds=lat_bnds)

ds = xr.open_mfdataset(sss_files, preprocess=partial_func)[["sos", "sos_error"]]
salt = ds.sos
ds.close()

ds = xr.open_mfdataset(sst_files, preprocess=partial_func)[["analysed_sst", "analysis_error"]]
temp = ds.analysed_sst - 273.15
ds.close()

def _preprocess(x, lon_bnds, lat_bnds):
    return x.sel(longitude=slice(*lon_bnds), latitude=slice(*lat_bnds))

partial_func = partial(_preprocess, lon_bnds=lon_bnds, lat_bnds=lat_bnds)

ds = xr.open_mfdataset(sla_files, preprocess=partial_func)[["adt"]]
zeta = ds.adt
ds.close()

print("Data loaded successfully.")


salt_mean = salt.mean(dim=('time', 'depth'))
zeta_mean = zeta.mean(dim='time')
temp_mean = temp.mean(dim='time')

salt = salt.sel(time=date).mean(dim='depth')
zeta = zeta.sel(time=date)
temp = temp.sel(time=date)


salt = salt.mean(dim='time')
temp = temp.mean(dim='time')

print("Data processed successfully.")

fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"MEASURES {date}")

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
label = 'SSS [psu]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
cb.set_ticks(np.arange(norm.vmin, norm.vmax + 0.5, 0.5))

# Free surface
ax = axs[1]
ax.set_title("Free Surface")

norm = mpl.colors.Normalize(vmin=0, vmax=1.5)
data = zeta
cmap = ssh_cmap
label = 'SSH [m]'
pcm = ax.pcolormesh(data.longitude, data.latitude, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Temperature
ax = axs[2]
ax.set_title("Temperature")

norm = mpl.colors.Normalize(vmin=20, vmax=30)
data = temp
cmap = sst_cmap
label = 'SST [°C]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')


fig.tight_layout()
filename = os.path.join(figures, f"physical_measures_{date}.png")
fig.savefig(filename, dpi=300)
print(f'{filename} saved.')
plt.show()



fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"MEASURES")

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
data = salt_mean
cmap = sss_cmap
label = r'$\overline{\overline{SSS}}$ [psu]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
cb.set_ticks(np.arange(norm.vmin, norm.vmax + 0.5, 0.5))

# Free surface
ax = axs[1]
ax.set_title("Free Surface")

norm = mpl.colors.Normalize(vmin=0, vmax=1.5)
data = zeta_mean
cmap = ssh_cmap
label = r'$\overline{\overline{SSH}}$ [m]'
pcm = ax.pcolormesh(data.longitude, data.latitude, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Temperature
ax = axs[2]
ax.set_title("Temperature")

norm = mpl.colors.Normalize(vmin=20, vmax=30)
data = temp_mean
cmap = sst_cmap
label = r'$\overline{\overline{SST}}$ [°C]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')


fig.tight_layout()
filename = os.path.join(figures, f'physical_measures_{year}_mean.png')
fig.savefig(filename, dpi=300)
print(f'{filename} saved.')
plt.show()



fig, axs = plt.subplots(1, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f"MEASURES {date}")

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

norm = mpl.colors.Normalize(vmin=-1, vmax=1)
data = salt - salt_mean
cmap = sss_cmap
label = r'SSA [psu]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmocean.cm.tarn, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')
cb.set_ticks(np.arange(norm.vmin, norm.vmax + 0.5, 0.5))

# Free surface
ax = axs[1]
ax.set_title("Free Surface")

norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
data = zeta - zeta_mean
cmap = ssh_cmap
label = r'SLA [m]'
pcm = ax.pcolormesh(data.longitude, data.latitude, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')

# Temperature
ax = axs[2]
ax.set_title("Temperature")

norm = mpl.colors.Normalize(vmin=-4, vmax=4)
data = temp - temp_mean
cmap = sst_cmap
label = r'STA [°C]'

pcm = ax.pcolormesh(data.lon, data.lat, data, cmap=cmocean.cm.balance, norm=norm, transform=ccrs.PlateCarree())
cb = plt.colorbar(pcm, ax=ax, label=label, orientation='vertical')


fig.tight_layout()
filename = os.path.join(figures, f"anomalies_measures_{date}.png")
fig.savefig(filename, dpi=300)
print(f'{filename} saved.')
plt.show()
