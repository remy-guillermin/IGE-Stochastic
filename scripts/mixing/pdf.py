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

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_avg_suf.nc'
#data_path = '/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_suf.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
#grid_path = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
glorys_path = '/home/guilremy/IGE-Stochastic/Data/glorys_avg_suf.nc'
#glorys_path = '/lus/work/CT1/c1601279/rguillermin/GLORYS/glorys_avg_suf.nc'
data_path2 = '/home/guilremy/IGE-Stochastic/Data/swio_his.nc'
dataset_path = '/home/guilremy/IGE-Stochastic/Data/datasets'

visc_cmap = cmocean.cm.amp
diff_cmap = cmocean.cm.tempo

ds = xr.open_dataset(data_path2)
AKv = ds.AKv
AKt = ds.AKt
AKs = ds.AKs
s_rho = ds.s_rho
Cs_rho = ds.Cs_rho
hc = ds.hc.values
ds.close()

date_index = -1
date = np.datetime_as_string(AKv.time.data[date_index], 'h')

s_rho_index = -1

fig, axes = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})

fig.suptitle(rf'Vertical Coefficients [m²/s] - {date} - $\sigma=${s_rho[s_rho_index]:.2f}')

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

plt.tight_layout()

plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/coefficients_{date}.png', dpi=300, transparent=True)
plt.show()



dates_index = (-36, -1)
dates = np.datetime_as_string(AKv.time.data[dates_index[0]: dates_index[1]], 'D')

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
fig.suptitle(rf'Vertical Coefficients PDF - {dates[0]} to {dates[-1]} - $\sigma=${s_rho[s_rho_index]:.2f}')

# Temperature Diffusion ($AKt$)
ax = axes[0]
ax.set_title('Temperature Diffusion ($AKt$)')

data = AKt[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
ax.hist(data, bins, color='blue', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKt$')

# Viscosity ($AKv$)
ax = axes[1]
ax.set_title('Viscosity ($AKv$)')

data = AKv[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
ax.hist(data, bins, color='green', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKv$')

# Salinity Diffusion ($AKs$)
ax = axes[2]
ax.set_title('Salinity Diffusion ($AKs$)')

data = AKs[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
n, _, _ = ax.hist(data, bins, color='red', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKs$')


for ax in axes:
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.legend()

axes[0].set_ylabel('Count')
plt.tight_layout()

plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/coefficients_PDF_{len(dates)}days.png', dpi=300, transparent=True)
plt.show()
