import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle

# Load data (same as your existing code)
data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

g = xr.open_dataset(grid)
h = g['h'][:, :]
lon = g['lon_rho'][:, :]  # Longitude
lat = g['lat_rho'][:, :]  # Latitude
angle = g['angle'][:, :]  # Deformation
msk = g['mask_rho'][:, :]  # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

# Define zones
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)]
names = ['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
colors = ['saddlebrown', 'darkorchid', 'navy', 'teal']

# Create figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Bathymetry plot
cmap = cmocean.cm.gray
a, b = 0, 5000
c = int((b-a)*2/1000+1)
levels = np.linspace(a, b, c*2-1)
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, h, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
cs = ax.contour(lon, lat, h, levels=np.arange(0, 6000, 1000), colors='k', linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs, fmt='%d', inline=True, fontsize=5)
cs2 = ax.contour(lon, lat, h, levels=np.arange(-300, 600, 400), colors='red', linestyles='dashed', linewidths=0.3, transform=ccrs.PlateCarree())
plt.clabel(cs2, fmt='%d', inline=True, fontsize=5)
ax.contour(lon, lat, msk, colors='k', linewidths=0.1)
ax.contourf(lon, lat, msk_inv, colors='lightgray')

# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 6, 'color': 'k'}
gl.ylabel_style = {'size': 6, 'color': 'k'}

# Add rectangles and labels
for (xmin, xmax, ymin, ymax), name, color in zip(boxes, names, colors):
    ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                           linewidth=2, edgecolor=color, facecolor='none', transform=ccrs.PlateCarree()))
    ax.text((xmin + xmax) / 2, ymax + 0.5, name, color=color, fontsize=8, ha='center', va='bottom', transform=ccrs.PlateCarree())

# Add bathymetry colorbar
cb = fig.colorbar(pcm, ax=ax, label='Bathymetry (m)')
colorbaryticks = np.linspace(a, b, c)
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=6)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=6)
text.set_font_properties(font)

plt.show()