# Libs
import cmcrameri
import cmocean
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from shapely.ops import unary_union
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import croco_plot as cplot
from matplotlib.patches import Rectangle


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

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_his.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'

lon, lat, _, _, msk, msk_inv, _, h = cplot.utils.load_grid(grid_path)

ds = xr.open_dataset(data_path)
AKv = ds.AKv
AKt = ds.AKt
AKs = ds.AKs
s_rho = ds.s_rho
Cs_rho = ds.Cs_rho
hc = ds.hc.values
ds.close()

depth_sigma = cplot.utils.calc_depth(s_rho, Cs_rho, hc, h)

date_index = -1
fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

date = np.datetime_as_string(AKv.time.data[date_index], 'h')

fig.suptitle('Vertical Temperature Diffusion Coefficient')

for ax, index in zip(axes, lon_index):
    ax.plot(lat[:, index].values, -h[:, index].values, marker='o', linestyle='-', color='k', markersize=1)
    ax.fill_between(lat[:, index], -h[:, index].values, y2=min(-h[:, index]), color='lightgrey')
    ax.set_title(f'Along {np.round(lon[0, index].values,2)} °E - {date}')
    data = AKv[date_index, :-1, :, index]
    data = data.where((data != 0), np.nan)
    
    offset = (ax.get_ylim()[1] - ax.get_ylim()[0])/20
    
    norm = mpl.colors.LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data))
    pcm = ax.pcolormesh(np.tile(lat[:, index].values, (50, 1)), depth_sigma[:, :, index], data, cmap=diff_cmap, norm=norm)
    cb = plt.colorbar(pcm, ax=ax, orientation='vertical', label='[m²/s]')
    
    for box, name, color in zip(boxes, names, colors):
        if lon[0, index].values < box[1] and lon[0, index].values > box[0]:
            ax.axvspan(box[2], box[3], color=color, alpha=0.05)
            ax.axvline(box[2], color=color, linewidth=2)
            ax.axvline(box[3], color=color, linewidth=2)
            ax.text(box[2] + (box[3] - box[2])/2, ax.get_ylim()[0] + offset, name, fontsize=11, verticalalignment='bottom', horizontalalignment='center')


for ax in axes:
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{tick:.0f}°N' for tick in ticks])
    
    ticks = ax.get_yticks()[1:-1]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick/1000)} km' for tick in ticks])
    
axes[-1].set_xlim((-35,6))

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/coefficients_lon_slices.png', dpi=300, transparent=True)

date_index = -1
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

date = np.datetime_as_string(AKv.time.data[date_index], 'h')

fig.suptitle('Vertical Temperature Diffusion Coefficient')

for ax, index in zip(axes, lat_index):
    ax.plot(lon[index, :].values, -h[index, :].values, marker='o', linestyle='-', color='k', markersize=1)
    ax.fill_between(lon[index, :], -h[index, :].values, y2=min(-h[index, :]), color='lightgrey')
    if lat[index, 0].values > 0:
        ax.set_title(f'Along {np.round(lat[index, 0].values,2)} °N - {date}')
    else:
        ax.set_title(f'Along {np.round(-lat[index, 0].values,2)} °S - {date}')
    data = AKv[date_index, :-1, index, :]
    data = data.where((data != 0), np.nan)
    
    norm = mpl.colors.LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data))
    pcm = ax.pcolormesh(np.tile(lon[index, :].values, (50, 1)), depth_sigma[:, index, :], data, cmap=diff_cmap, norm=norm)
    cb = plt.colorbar(pcm, ax=ax, orientation='vertical', label='[m²/s]')
    
    offset = (ax.get_ylim()[1] - ax.get_ylim()[0])/20
    
    for box, name, color in zip(boxes, names, colors):
        if lat[index, 0].values < box[3] and lat[index, 0].values > box[2]:
            ax.axvspan(box[0], box[1], color=color, alpha=0.05)
            ax.axvline(box[0], color=color, linewidth=2)
            ax.axvline(box[1], color=color, linewidth=2)
            if name == 'Mayotte-Comores':
                ax.text(box[0] + (box[1] - box[0])/2, ax.get_ylim()[0] + offset, name, fontsize=10, verticalalignment='bottom', horizontalalignment='center')
            else:
                ax.text(box[0] + (box[1] - box[0])/2, ax.get_ylim()[1] - offset, name, fontsize=10, verticalalignment='top', horizontalalignment='center')

for ax in axes:
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{tick:.0f}°E' for tick in ticks])
    
    ticks = ax.get_yticks()[1:-1]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{-int(tick/1000)} km' for tick in ticks])

 
axes[-1].set_xlim((30, 65))

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/coefficients_lat_slices.png', dpi=300, transparent=True)
plt.show()


