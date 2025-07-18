# LIBRAIRIES
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import cmcrameri
import cartopy.crs as ccrs
import glob
import os
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# PATH
figures = '/lus/home/CT1/c1601279/rguillermin/IGE-Stochastic/figures'
store = '/lus/store/CT1/c1601279/lweiss/run_croco/SWIOSE_dev/'
stoens = 'run_swio2_stoens30_2018_CD'
grid = '/lus/store/CT1/c1601279/lweiss/grid/croco_grid_swio2.nc'

SWIO = (25, 69, -36, 7)
vort_cmap=cmcrameri.cm.vik
velo_cmap=cmocean.cm.speed
ener_cmap=cmcrameri.cm.devon
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

# Define zones
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -30, -21), (52, 60, -24, -16)]
names = ['Equator', 'Islands', 'South-Moz', 'Mascarene']
cmap = cmocean.cm.tempo_r
n_colors = 4
box_colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
box_colors[0] = mcolors.to_rgba('palevioletred')
box_colors[-1] = mcolors.to_rgba('sandybrown')

lat_index = [430, 280]
points = [(280, 380), (280, 190), (430, 180)]
points_names = ['Mascarene', 'Islands', 'Somali']

def calc_depth(
    s: xr.DataArray, 
    Cs: xr.DataArray, 
    hc: float, 
    h: xr.DataArray
) -> xr.DataArray:
    """
    Compute depth using the S-coordinate transformation.

    Parameters
    ----------
    s : xr.DataArray
        S-coordinate at RHO-points, typically ranging from -1 (surface) to 0 (bottom).
    Cs : xr.DataArray
        S-coordinate stretching curves at W-points, defining the vertical stretching.
    hc : float
        Critical depth parameter (in meters), influencing vertical terrain-following transformation.
    h : xr.DataArray
        Bathymetric depth at RHO-points (in meters), representing the seafloor depth.

    Returns
    -------
    xr.DataArray
        Computed depth at RHO-points.
    """
    N = len(s)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
        depth[k, :, :] = z0[k, :, :] * h
    return depth

g = xr.open_dataset(grid)[['lon_rho', 'lat_rho', 'mask_rho', 'h']]
lon = g.lon_rho
lat = g.lat_rho
eta_rho = g.eta_rho
xi_rho = g.xi_rho
h = g.h
g.close()

ensemble_str_files = [os.path.join(store, stoens, f"{i:03d}swiose_avg.nc") for i in range(1, 31)]

ds = xr.open_dataset(ensemble_str_files[0])
depth_sigma = calc_depth(ds.s_rho, ds.Cs_rho, ds.hc, h)
s_rho = ds.s_rho
temp = ds.temp
ds.close()

depth_level = xr.Dataset(
    data_vars=dict(
        depth_level=(['s_rho', 'eta_rho', 'xi_rho'], depth_sigma),
    ),
    coords=dict(
        eta_rho=eta_rho,
        xi_rho=xi_rho,
        s_rho=s_rho,
    )
)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.set_extent(SWIO)   
ax.coastlines(resolution='50m')
ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

# Bathymetry plot
cmap = cmocean.cm.gray
a, b = 0, 5000
c = int((b-a)*2/1000+1)
levels = np.linspace(a, b, c*2-1)
norm = mcolors.BoundaryNorm(levels, cmap.N)

pcm = ax.pcolormesh(lon, lat, h, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())


# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='-', linewidth=0.2, color='k', alpha=0.5)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'k'}
gl.ylabel_style = {'size': 10, 'color': 'k'}

# Add rectangles and labels
# for (xmin, xmax, ymin, ymax), name, color in zip(boxes, names, box_colors):
#     if name == 'Equator'or name == 'Islands':
#         ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
#                             linewidth=2, edgecolor=color, facecolor='none', transform=ccrs.PlateCarree()))
#         ax.text((xmin + xmax) / 2, ymax + 0.5, name, color=color, fontsize=10, ha='center', va='bottom', transform=ccrs.PlateCarree())

# i = 1
# for index in lat_index:
#     ax.plot(lon[index, :], lat[index, :], 'r-',)
#     ax.text(lon[index, -1] + 0.5, lat[index, -1], f'S{i}', color='red', ha='left', va='center')
#     i += 1


for k, (i, j) in enumerate(points):
    ax.scatter(lon[i, j], lat[i, j], 10, color='red')
    ax.text(lon[i, j], lat[i, j] - 0.5, points_names[k], color='red', ha='center', va='top')

fig.tight_layout()
fig.savefig(f'{figures}/bathy_points.png', dpi=300)
plt.show()
