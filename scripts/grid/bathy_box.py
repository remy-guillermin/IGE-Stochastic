import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle

# Load data (same as your existing code)
# LOCAL
# grid = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
# figures = '/home/guilremy/IGE-Stochastic/figures'

# CLUSTER
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'

g = xr.open_dataset(grid)
h = g['h'][:, :] # Bathymetry
lon = g['lon_rho'][:, :]  # Longitude
lat = g['lat_rho'][:, :]  # Latitude
angle = g['angle'][:, :]  # Deformation
msk = g['mask_rho'][:, :]  # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()

# Define zones
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
colors = ['saddlebrown', 'darkorchid', 'navy', 'teal']

# Create figure
fig = plt.figure(figsize=(10, 8))
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

ax.add_feature(minor_islands)

# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k', zorder=4)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'color': 'k'}
gl.ylabel_style = {'color': 'k'}

# Add rectangles and labels
for (xmin, xmax, ymin, ymax), name, color in zip(boxes, names, colors):
    ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                           linewidth=2, edgecolor=color, facecolor='none', transform=ccrs.PlateCarree()))
    ax.text((xmin + xmax) / 2, ymax + 0.5, name, color=color, fontsize=10, ha='center', va='bottom', transform=ccrs.PlateCarree())

# Add bathymetry colorbar
cb = fig.colorbar(pcm, ax=ax, label='Bathymetry (m)')
colorbaryticks = np.linspace(a, b, c)
posax = ax.get_position()
poscb = cb.ax.get_position()
cb.ax.set_position([0.76, posax.y0, poscb.width, posax.height])
cb.set_ticks(colorbaryticks)
cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=10)
text = cb.ax.yaxis.label
font = mpl.font_manager.FontProperties(size=10)
text.set_font_properties(font)

fig.tight_layout()
fig.savefig(f'{figures}/bathy_zones.png', dpi=300, transparent=True)
plt.show()
