import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cartopy.crs as ccrs
import cmcrameri
from scipy.ndimage import label as ndi_label, sum as ndi_sum
from skimage.measure import regionprops, find_contours
import croco_plot as cplot
#%matplotlib tk

#---------------------------------#
# PARAMETERS AND DATA LOADING
#---------------------------------#

output_dir = 'Frames'
os.makedirs(output_dir, exist_ok=True)

# Constants
vort_cmap = cmcrameri.cm.vik
figsize = (14, 12)
min_vort = 1.5e-5
min_area = 50
min_circularity = 0.25
gridline_style = {'draw_labels': False, 'linestyle': '--', 'linewidth': 0.3}

# Load data
data_path = '/home/guilremy/IGE-Stochastic/Data/swio_avg_suf.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
lon, lat, pm, pn, msk, msk_inv, angle, _ = cplot.utils.load_grid(grid_path)
velocity, vorticity = cplot.utils.load_data(data_path, ('V_hor', 'vort'))
velocity, vorticity = velocity[:, 0, :, :], vorticity[:, 0, :, :]

#---------------------------------#
# PLOTTING SETUP (Only once!)
#---------------------------------#

fig, axes = plt.subplots(2, 2, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
plt.subplots_adjust(wspace=0.05, hspace=0.15)

mask_cmap = mcolors.ListedColormap(["white", "red"])
mask_norm = mcolors.BoundaryNorm([0, 0.5, 1], mask_cmap.N)
levels = np.linspace(-0.00005, 0.00005, 51)
norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)

# --- Panels ---
pcm = axes[0, 0].pcolormesh(lon, lat, np.zeros_like(vorticity[0]), cmap=vort_cmap, norm=norm, transform=ccrs.PlateCarree())
pcm_mask = axes[0, 1].pcolormesh(lon, lat, np.zeros_like(vorticity[0], dtype=bool), cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())
pcm_size = axes[1, 0].pcolormesh(lon, lat, np.zeros_like(vorticity[0], dtype=bool), cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())
pcm_shape = axes[1, 1].pcolormesh(lon, lat, np.zeros_like(vorticity[0], dtype=bool), cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

# Shared coastlines, contours, gridlines
for a in axes.flat:
    a.coastlines(zorder=100)
    a.contourf(lon, lat, msk_inv, colors='lightgray', transform=ccrs.PlateCarree())
    a.contour(lon, lat, msk, colors='k', linewidths=0.2, transform=ccrs.PlateCarree())
    gl = a.gridlines(crs=ccrs.PlateCarree(), **gridline_style)
    gl.bottom_labels = False
    gl.left_labels = False

# Colorbars (set only once)
fig.colorbar(pcm, ax=axes[0, 0], orientation='vertical', fraction=0.046, pad=0.04, label='Vorticity [m²/s]')
for ax, label in zip(axes.flat[1:], [f'> {min_vort}', f'> {min_vort}, A>{min_area}', f'> {min_vort}, A>{min_area}, C>{min_circularity}']):
    cb = fig.colorbar(pcm_mask, ax=ax, orientation='vertical', fraction=0.046, pad=0.04, ticks=[0, 1], label=label)
    cb.ax.set_yticklabels(['No', 'Yes'])

date = np.datetime_as_string(vorticity.time.data[0], 'D')
fig.suptitle(f"Vorticity Analysis and Filtering | Date: {date}", fontsize=14)
plt.tight_layout()
   
#---------------------------------#
# MAIN LOOP FOR TIME STEPS
#---------------------------------#

for idx in range(vorticity.shape[0]):
    vort_values = vorticity[idx].data
    mask_values = vort_values > min_vort
    labeled_mask, num_features = ndi_label(mask_values)

    # --- Size filtering ---
    region_areas = ndi_sum(mask_values, labeled_mask, index=range(1, num_features + 1))
    valid_labels_size = np.where(region_areas >= min_area)[0] + 1
    filtered_mask_size = np.isin(labeled_mask, valid_labels_size)

    # --- Shape (circularity) filtering ---
    filtered_mask_shape = np.zeros_like(mask_values, dtype=bool)
    for prop in regionprops(labeled_mask):
        area, perimeter = prop.area, prop.perimeter
        if area >= min_area and (4 * np.pi * area) / (perimeter ** 2) >= min_circularity:
            filtered_mask_shape[labeled_mask == prop.label] = True

    # --- Panel Updates ---
    pcm.set_array(vort_values.ravel())
    pcm_mask.set_array(mask_values.ravel())
    pcm_size.set_array(filtered_mask_size.ravel())
    pcm_shape.set_array(filtered_mask_shape.ravel())

    # --- Titles and figure title ---
    date = np.datetime_as_string(vorticity.time.data[idx], 'D')
    axes[0, 0].set_title("Vorticity")
    axes[0, 1].set_title(rf"Value filter ($|\omega|>{min_vort}$) | $N_f={num_features}$")
    axes[1, 0].set_title(rf"Size filter ($A>{min_area}$) | $N_f={len(valid_labels_size)}$")
    axes[1, 1].set_title(rf"Shape filter ($A>{min_area}$, $C>{min_circularity}$) | $N_f={np.unique(ndi_label(filtered_mask_shape)[0]).size - 1}$")
    fig.suptitle(f"Vorticity Analysis and Filtering | Date: {date}", fontsize=14)

    # --- Save Frame ---
    output_filename = os.path.join(output_dir, f'frame_{idx:04d}_{date}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f'Saved {output_filename}')

for prop in regionprops(labeled_mask):
    area, perimeter = prop.area, prop.perimeter
    circularity = (4 * np.pi * area) / (perimeter ** 2) if perimeter > 0 else 0
    y, x = map(int, prop.centroid)
    Lon, Lat = lon.data[y, x], lat.data[y, x]
    
    if area >= min_area:
        # Print result
        print(
            f'Perimeter: {perimeter:.2f}, Area: {area:.0f}, Circularity: {circularity:.2f}, '
            f'Centroid: ({Lon:.2f}, {Lat:.2f}), '
        )
