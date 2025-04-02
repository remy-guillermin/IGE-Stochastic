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

fname = '/home/guilremy/Téléchargements/gadm41_IND_shp/gadm41_IND_1.shp'

adm1_shapes = list(shpreader.Reader(fname).geometries())

SWI = (72, 86, 5, 20)

print(len(adm1_shapes))

shapes = [shape for shape in adm1_shapes if shape.bounds[3] < SWI[3]]

print(len(shapes))

fname = '/home/guilremy/Téléchargements/gadm41_IND_shp/gadm41_IND_2.shp'

adm2_shapes = list(shpreader.Reader(fname).geometries())

print(len(adm2_shapes))

shapes2 = [shape for shape in adm2_shapes if shape.bounds[3] < SWI[3]]

print(len(shapes2))

LAK = shapes[5]
KER = shapes[4]
TN = shapes[7]
TEL = shapes[8]

msk = shapes[1].contains(shapes2)
CAPRYS = [shape for shape in shapes2 if msk[shapes2.index(shape)]]
print(len(CAPRYS))

msk = shapes[3].contains(shapes2)
SIKNIKCK = [shape for shape in shapes2 if msk[shapes2.index(shape)]]
print(len(SIKNIKCK))

SIK = unary_union([SIKNIKCK[1], SIKNIKCK[2], SIKNIKCK[4], SIKNIKCK[7], SIKNIKCK[8], SIKNIKCK[9], SIKNIKCK[10], SIKNIKCK[12], SIKNIKCK[16], SIKNIKCK[18], SIKNIKCK[19], SIKNIKCK[21], SIKNIKCK[22], SIKNIKCK[24], SIKNIKCK[25], SIKNIKCK[26]])
NIK = unary_union([SIKNIKCK[0], SIKNIKCK[3], SIKNIKCK[5], SIKNIKCK[6], SIKNIKCK[13], SIKNIKCK[14], SIKNIKCK[15], SIKNIKCK[17], SIKNIKCK[20], SIKNIKCK[23], SIKNIKCK[29]])
CK = unary_union([shapes[2], SIKNIKCK[11], SIKNIKCK[27], SIKNIKCK[28]])

RYS = unary_union([CAPRYS[0], CAPRYS[1], CAPRYS[5], CAPRYS[12]])
CAP = unary_union([shape for shape in CAPRYS if not RYS.contains(shape)])

fig, ax = plt.subplots(1, 1, figsize=(6, 6), subplot_kw={'projection': ccrs.PlateCarree()})

fig.suptitle('Sub division of the South West India')

ax.set_extent(SWI)   
ax.coastlines(resolution='10m')
ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black', zorder=3)
ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5, zorder=3)
ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5, zorder=3)

ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='black', facecolor='none', zorder=4)

ax.add_geometries(LAK.envelope, ccrs.PlateCarree(), edgecolor='black', facecolor='teal', zorder=5)
ax.add_geometries(RYS, ccrs.PlateCarree(), edgecolor='black', facecolor='green', zorder=5)
ax.add_geometries(CAP, ccrs.PlateCarree(), edgecolor='black', facecolor='green', zorder=5)
ax.add_geometries(KER, ccrs.PlateCarree(), edgecolor='black', facecolor='green', zorder=5)
ax.add_geometries(TN, ccrs.PlateCarree(), edgecolor='black', facecolor='green', zorder=5)
ax.add_geometries(TEL, ccrs.PlateCarree(), edgecolor='black', facecolor='green', zorder=5)
ax.add_geometries(NIK, ccrs.PlateCarree(), edgecolor='black', facecolor='teal', zorder=5)
ax.add_geometries(SIK, ccrs.PlateCarree(), edgecolor='black', facecolor='teal', zorder=5)
ax.add_geometries(CK, ccrs.PlateCarree(), edgecolor='black', facecolor='teal', zorder=5)

ax.text(LAK.centroid.x, LAK.centroid.y, 'LAK', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(TEL.centroid.x, TEL.centroid.y, 'TEL', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(SIK.centroid.x, SIK.centroid.y, 'SIK', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(TN.centroid.x, TN.centroid.y, 'TN', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(NIK.centroid.x, NIK.centroid.y, 'NIK', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(CAP.centroid.x, CAP.centroid.y, 'CAP', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(RYS.centroid.x, RYS.centroid.y, 'RYS', fontsize=10, verticalalignment='center', horizontalalignment='center', zorder=6)
ax.text(CK.centroid.x-0.2, CK.centroid.y+0.5, 'CK', fontsize=10, verticalalignment='center', horizontalalignment='center', rotation=-60, zorder=6)
ax.text(KER.centroid.x, KER.centroid.y, 'KER', fontsize=10, verticalalignment='center', horizontalalignment='center', rotation=-60,zorder=6)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.2, color='k')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'color': 'k'}
gl.ylabel_style = {'color': 'k'}

plt.tight_layout()
plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/India_subdiv.png', dpi=300, transparent=True)
plt.show()
