# Libs
import cmcrameri
import cmocean
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.ops import unary_union
import numpy as np
import xarray as xr
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

grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
data_path = '/home/guilremy/IGE-Stochastic/Data/swio_his.nc'

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

# Create figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

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

# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'k'}
gl.ylabel_style = {'size': 10, 'color': 'k'}

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

for index in lon_index:
    ax.plot(lon[:, index], lat[:, index], 'r-')
    
for index in lat_index:
    ax.plot(lon[index, :], lat[index, :], 'r-')

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/bathy_slices.png', dpi=300, transparent=True)
plt.show()
