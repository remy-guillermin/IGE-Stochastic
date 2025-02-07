import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cmocean
import numpy as np
import xarray as xr
from PIL import Image

scratch = '/lus/scratch/CT1/c1601279/lweiss/CROCO/SWIOSE/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swiose_0002'

# Charger le fichier netcdf
ds = xr.open_dataset(scratch + simu + '/swiose_his.nc')

temp = ds['temp'][:, :, :, :]  # Sélectionne le premier niveau s_rho
lon = ds['lon_rho'][:,:]
lat = ds['lat_rho'][:,:]
msk = ds['mask_rho'][:,:]
ds.close()

# Inverser le masque pour le rendre transparent où msk==0
msk_inv = np.where(msk==0, msk, np.nan)

# Coordonnées des points à zoomer
# point1 = (-1, -1, 249, 195)
# point2 = (-1, 0, 249, 195)
point1 = (-1, -1, 302, 90)
point2 = (-1, 0, 302, 90)

S = 40

# Définition des indices de grille pour le zoom
zoom1 = (slice(point1[2]-int(S/2), point1[2]+int(S/2)), slice(point1[3]-int(S/2), point1[3]+int(S/2)))
zoom2 = (slice(point2[2]-int(S/2), point2[2]+int(S/2)), slice(point2[3]-int(S/2), point2[3]+int(S/2)))

# Création des figures
cmap = cmocean.cm.thermal
a = 10
b = 31
c = int((b-a)+1)
levels = np.linspace(a, b, c*2-1)
norm = mpl.colors.BoundaryNorm(levels, cmap.N)

fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})

for i, (point, zoom, ax) in enumerate(zip([point1, point2], [zoom1, zoom2], axes)):

    pcm = ax.pcolormesh(lon[zoom], lat[zoom], temp[point[0], point[1], zoom[0], zoom[1]], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.contour(lon[zoom], lat[zoom], msk[zoom], colors='k', linewidths=0.1)
    ax.contourf(lon[zoom], lat[zoom], msk_inv[zoom], colors='lightgray')
    ax.set_title(f'Zoom autour de {simu} point {point}', size=8)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.4)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8, 'color': 'k'}
    gl.ylabel_style = {'size': 8, 'color': 'k'}

    cb = fig.colorbar(pcm, ax=ax, label='Sea Surface Temperature °C')
    colorbaryticks = np.linspace(a, b, c)
    posax = ax.get_position()
    poscb = cb.ax.get_position()
    # cb.ax.set_position([0.92 + i*0.02, posax.y0, poscb.width, posax.height])
    cb.ax.set_position([0.3 + posax.x0, posax.y0, poscb.width, posax.height])
    cb.set_ticks(colorbaryticks)
    cb.ax.set_yticklabels(colorbaryticks.astype(int), fontsize=8)
    text = cb.ax.yaxis.label
    font = mpl.font_manager.FontProperties(size=8)
    text.set_font_properties(font)

# Sauvegarde des figures
plt.savefig(path_fig + 'zoomed_maps_max_' + str(point[2]) + '_' + str(point[3]) + simu + '.png', dpi=300, bbox_inches='tight')

