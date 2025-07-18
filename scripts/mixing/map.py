# Libs
import cmcrameri
import cmocean
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import os

SWIO = (25, 69, -36, 7)

vort_cmap=cmcrameri.cm.vik
velo_cmap=cmocean.cm.speed
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

# Define zones
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
colors = ['saddlebrown', 'darkorchid', 'navy', 'teal']

lat_index = (172, 280, 430)
lon_index = (160, 335)

visc_cmap = cmocean.cm.amp
diff_cmap = cmocean.cm.tempo

# LOCAL
# grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
# figures = '/home/guilremy/IGE-Stochastic/figures'

# CLUSTER
grid_path = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
data_path = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/data/swio_his.nc'

g = xr.open_dataset(grid_path)
h = g['h'][:, :]
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

ds = xr.open_dataset(data_path)
AKv = ds.AKv
AKt = ds.AKt
AKs = ds.AKs
s_rho = ds.s_rho
Cs_rho = ds.Cs_rho
hc = ds.hc.values
ds.close()

fill_value = 9.96921e+36
AKv = AKv.where((AKv != fill_value), np.nan)
AKt = AKt.where((AKt != fill_value), np.nan)
AKs = AKs.where((AKs != fill_value), np.nan)

date_index = -1
date = np.datetime_as_string(AKv.isel(time=date_index).time, 'h')

s_rho_index = -1

fig, axes = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})

fig.suptitle(rf'Vertical Coefficients [m²/s] - {date} - $\sigma=${s_rho[s_rho_index]:.2f}')

ax = axes[0]
ax.set_title('Temperature Diffusion ($AKt$)')

data = AKt[date_index,s_rho_index,:,:]
data = data.where((data != 0), np.nan)
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

norm = mpl.colors.LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data))

pcm_t = ax.pcolormesh(lon[:,:], lat[:,:], data, cmap=diff_cmap, norm=norm, transform=ccrs.PlateCarree())
cb_t = plt.colorbar(pcm_t, ax=ax, orientation='vertical', label='[m²/s]')

ax = axes[1]
ax.set_title('Viscosity ($AKv$)')

data = AKv[date_index,s_rho_index,:,:]
data = data.where((data != 0), np.nan)
print(f'AKv: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

norm = mpl.colors.LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data))

pcm_v = ax.pcolormesh(lon[:,:], lat[:,:], data, cmap=visc_cmap, norm=norm, transform=ccrs.PlateCarree())
cb_v = plt.colorbar(pcm_v, ax=ax, orientation='vertical', label='[m²/s]')

ax = axes[2]
ax.set_title('Salinity Diffusion ($AKs$)')

data = AKs[date_index,s_rho_index,:,:]
data = data.where((data != 0), np.nan)
print(f'AKs: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

norm = mpl.colors.LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data))
pcm_v = ax.pcolormesh(lon[:,:], lat[:,:], data, cmap=diff_cmap, norm=norm, transform=ccrs.PlateCarree())
cb_v = plt.colorbar(pcm_v, ax=ax, orientation='vertical', label='[m²/s]')

for ax in axes:
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
    gl.xlabel_style = {'color': 'k'}
    gl.ylabel_style = {'color': 'k'}
    

fig.tight_layout()

filename = f'{figures}/coefficients_{date}.png'
fig.savefig(filename, dpi=300)
print(f'Saved in {filename}')
plt.show()
