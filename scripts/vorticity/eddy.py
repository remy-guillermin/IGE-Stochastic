import cmocean
import cmcrameri
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from scipy.ndimage import label, sum as ndi_sum
from skimage.measure import regionprops
import croco_plot as cplot

vort_cmap=cmcrameri.cm.vik

figsize=(14, 12)
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_avg_suf.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'

lon, lat, pm, pn, msk, msk_inv, angle, _ = cplot.utils.load_grid(grid_path)
velocity, vorticity = cplot.utils.load_data(data_path, ('V_hor','vort'))

velocity, vorticity = velocity[:,0,:,:], vorticity[:,0,:,:]

def update_plot(index):
	mask_values = np.abs(vorticity[index].data) > min_vort
	# Label connected components
	labeled_mask, num_features = label(mask_values)

	# Compute the area of each labeled region
	region_areas = ndi_sum(mask_values, labeled_mask, index=range(1, num_features + 1))

	# Create a new mask keeping only regions above the threshold
	filtered_mask = np.zeros_like(mask_values, dtype=bool)
	for i, area in enumerate(region_areas):
		if area >= min_area:
			filtered_mask[labeled_mask == (i + 1)] = True

	num_size_filtered = np.unique(labeled_mask[filtered_mask]).size
	if 0 in np.unique(labeled_mask[filtered_mask]):
		num_size_filtered -= 1

	# Compute properties of each labeled region
	filtered_mask2 = np.zeros_like(mask_values, dtype=bool)

	props = regionprops(labeled_mask)
	for prop in props:
		area = prop.area
		perimeter = prop.perimeter

		# Compute circularity (avoid division by zero)
		circularity = (4 * np.pi * area) / (perimeter**2) if perimeter > 0 else 0

		# Keep only large and round structures
		if area >= min_area and circularity >= min_circularity:
			filtered_mask2[labeled_mask == prop.label] = True

		num_shape_filtered = np.unique(labeled_mask[filtered_mask2]).size
		if 0 in np.unique(labeled_mask[filtered_mask2]):
			num_shape_filtered -= 1

	pcm.set_array(vorticity[index])
	pcmmask.set_array(mask_values)
	pcmfilt.set_array(filtered_mask)
	pcmfilt2.set_array(filtered_mask2)
	date = np.datetime_as_string(vorticity.time.data[index], 'D')
	ax.set_title(f"Vorticity SWIO {date}")
	ax2.set_title(rf"Vorticity masked SWIO {date} | $N_f={num_features}$")
	ax3.set_title(rf"Vorticity filtered (size) SWIO {date} | $N_f={num_size_filtered}$")
	ax4.set_title(rf"Vorticity filtered (shape) SWIO {date} | $N_f={num_shape_filtered}$")
	plt.draw()

date = np.datetime_as_string(vorticity.time.data[0], 'D')

fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_title(f"Vorticity SWIO {date}")

levels = np.linspace(-0.00005, 0.00005, 51)
norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)

pcm = ax.pcolormesh(lon[:, :], lat[:, :],  vorticity[0], cmap=vort_cmap, norm=norm, transform=ccrs.PlateCarree())

ax.contourf(lon, lat, msk_inv, colors='lightgray')
ax.contour(lon, lat, msk, colors='k', linewidths=0.2)

gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl.bottom_labels = False
gl.left_labels = False
gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}

cb = plt.colorbar(pcm, ax=ax, label=r'Vorticity [m²/s]', orientation='vertical')
cb.set_ticks(cb.get_ticks())
cb.ax.set_yticklabels(np.round(cb.get_ticks(), 6))

plt.tight_layout()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MASKING
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

min_vort = 1e-5

# Define custom colormap (white for False, red for True)
mask_cmap = mcolors.ListedColormap(["white", "red"])
mask_norm = mcolors.BoundaryNorm([0, 0.5, 1], mask_cmap.N)

# Boolean mask where condition is met
mask_values = np.abs(vorticity[0].data) > min_vort

labeled_mask, num_features = label(mask_values)

fig2, ax2 = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
ax2.set_title(rf"Vorticity masked SWIO {date} | $N_f={num_features}$")

pcmmask = ax2.pcolormesh(lon[:, :], lat[:, :],  mask_values, cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

ax2.contourf(lon, lat, msk_inv, colors='lightgray')
ax2.contour(lon, lat, msk, colors='k', linewidths=0.2)

gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl2.bottom_labels = False
gl2.left_labels = False
gl2.xlabel_style = gl2.ylabel_style = {'size': 8, 'color': 'k'}

cb2 = plt.colorbar(pcmmask, ax=ax2, label=r'$Vorticity > 1e-5$', orientation='vertical')
cb2.set_ticks(cb2.get_ticks())
cb2.ax.set_yticklabels(['No', '', 'Yes'])

plt.tight_layout()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIZE FILTERING
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

min_area = 50

# Label connected components
labeled_mask, num_features = label(mask_values)

# Compute the area of each labeled region
region_areas = ndi_sum(mask_values, labeled_mask, index=range(1, num_features + 1))

# Create a new mask keeping only regions above the threshold
filtered_mask = np.zeros_like(mask_values, dtype=bool)
for i, area in enumerate(region_areas):
    if area >= min_area:
        filtered_mask[labeled_mask == (i + 1)] = True

num_size_filtered = np.unique(labeled_mask[filtered_mask]).size
if 0 in np.unique(labeled_mask[filtered_mask]):
	num_size_filtered -= 1

fig3, ax3 = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
ax3.set_title(rf"Vorticity filtered (size) SWIO {date} | $N_f={num_size_filtered}$")

pcmfilt = ax3.pcolormesh(lon[:, :], lat[:, :],  filtered_mask, cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

ax3.contourf(lon, lat, msk_inv, colors='lightgray')
ax3.contour(lon, lat, msk, colors='k', linewidths=0.2)

gl3 = ax3.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl3.bottom_labels = False
gl3.left_labels = False
gl3.xlabel_style = gl3.ylabel_style = {'size': 8, 'color': 'k'}

cb3 = plt.colorbar(pcmfilt, ax=ax3, label=rf'$Vorticity > 1e-5$ & $N > {min_area}$', orientation='vertical')
cb3.set_ticks(cb3.get_ticks())
cb3.ax.set_yticklabels(['No', '', 'Yes'])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SHAPE FILTERING
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

min_circularity = 0.2

# Compute properties of each labeled region
filtered_mask = np.zeros_like(mask_values, dtype=bool)

props = regionprops(labeled_mask)
for prop in props:
    area = prop.area
    perimeter = prop.perimeter

    # Compute circularity (avoid division by zero)
    circularity = (4 * np.pi * area) / (perimeter**2) if perimeter > 0 else 0

    # Keep only large and round structures
    if area >= min_area and circularity >= min_circularity:
        filtered_mask[labeled_mask == prop.label] = True
        
num_shape_filtered = np.unique(labeled_mask[filtered_mask]).size
if 0 in np.unique(labeled_mask[filtered_mask]):
	num_shape_filtered -= 1

fig4, ax4 = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
ax4.set_title(rf"Vorticity filtered (shape) SWIO {date} | $N_f={num_shape_filtered}$")

pcmfilt2 = ax4.pcolormesh(lon[:, :], lat[:, :],  filtered_mask, cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

ax4.contourf(lon, lat, msk_inv, colors='lightgray')
ax4.contour(lon, lat, msk, colors='k', linewidths=0.2)

gl4 = ax4.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
gl4.bottom_labels = False
gl4.left_labels = False
gl4.xlabel_style = gl4.ylabel_style = {'size': 8, 'color': 'k'}

cb4 = plt.colorbar(pcmfilt2, ax=ax4, label=rf'$Vorticity > 1e-5$ & $N > {min_area}$ & $C > {min_circularity}$', orientation='vertical')
cb4.set_ticks(cb4.get_ticks())
cb4.ax.set_yticklabels(['No', '', 'Yes'])
