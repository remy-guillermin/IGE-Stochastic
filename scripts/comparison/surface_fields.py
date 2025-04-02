# Libs
import cmcrameri
import cmocean
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import croco_plot as cplot

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

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_avg_suf.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
glorys_path = '/home/guilremy/IGE-Stochastic/Data/glorys_avg_suf.nc'


vort_g, velo_g = cplot.utils.load_data(glorys_path, ('vort','V_hor'))
vort_c, velo_c = cplot.utils.load_data(data_path, ('vort','V_hor'))
lon, lat, pm, pn, msk, msk_inv, angle, h = cplot.utils.load_grid(grid_path)

mask = np.zeros_like(msk, dtype=bool)
mask[msk == 0] = True
mask[msk == 1] = False

velo_c, vort_c, velo_g, vort_g = velo_c[0, 0, :, :], vort_c[0, 0, :, :], velo_g[0, 0, :, :], vort_g[0, 0, :, :]
velo_c.data[mask] = np.nan
vort_c.data[mask] = np.nan

date = np.datetime_as_string(vort_g.time.data, 'D')

fig, axes = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={'projection': ccrs.PlateCarree()})

fig.suptitle('Surface velocity [m/s]')

for ax in axes:
    ax.set_extent(SWIO)   
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'color': 'k'}
    gl.ylabel_style = {'color': 'k'}

ax = axes[0]
ax.set_title('GLORYS')

pcm_g = ax.pcolormesh(velo_g.longitude, velo_g.latitude, velo_g, cmap=velo_cmap, norm=mpl.colors.LogNorm(vmin=0.1, vmax=2), transform=ccrs.PlateCarree())
cb_g = plt.colorbar(pcm_g, ax=ax, orientation='vertical')

ax = axes[1]
ax.set_title('CROCO')

pcm_c = ax.pcolormesh(lon.data, lat.data, velo_c, cmap=velo_cmap, norm=mpl.colors.LogNorm(vmin=0.1, vmax=2), transform=ccrs.PlateCarree())
cb_c = plt.colorbar(pcm_c, ax=ax, orientation='vertical')

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/velocity_dual_{date}.png', dpi=300, transparent=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={'projection': ccrs.PlateCarree()})

fig.suptitle(f'Surface vorticity [1/s] - {date}')

for ax in axes:
    ax.set_extent(SWIO)   
    ax.coastlines(resolution='50m')
    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'color': 'k'}
    gl.ylabel_style = {'color': 'k'}

ax = axes[0]
ax.set_title('GLORYS')

pcm_g = ax.pcolormesh(vort_g.longitude, vort_g.latitude, vort_g, cmap=vort_cmap, vmin=-0.0001, vmax=0.0001, transform=ccrs.PlateCarree())
cb_g = plt.colorbar(pcm_g, ax=ax, orientation='vertical')
ticks = cb_g.get_ticks()
cb_g.set_ticks(ticks)
cb_g.set_ticklabels([f'{tick:.1e}' for tick in ticks])

ax = axes[1]
ax.set_title('CROCO')

pcm_c = ax.pcolormesh(lon.data, lat.data, vort_c, cmap=vort_cmap, vmin=-0.0001, vmax=0.0001, transform=ccrs.PlateCarree())
cb_c = plt.colorbar(pcm_c, ax=ax, orientation='vertical')
ticks = cb_c.get_ticks()
cb_c.set_ticks(ticks)
cb_c.set_ticklabels([f'{tick:.1e}' for tick in ticks])

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/vorticity_dual_{date}.png', dpi=300, transparent=True)
plt.show()
