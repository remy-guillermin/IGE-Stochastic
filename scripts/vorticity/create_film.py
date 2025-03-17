# Libs
import cmcrameri
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
from scipy.ndimage import label, sum as ndi_sum
from skimage.measure import regionprops
import croco_plot as cplot

# Params

min_vort = 1.5e-5
min_area = 50
min_circularity = 0.2

SWIO = (25, 69, -36, 7)

vort_cmap=cmcrameri.cm.vik
gridline_style = {'draw_labels': True, 'linestyle': '--', 'linewidth': 0.3}

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_avg_suf.nc'
#data_path = '/lus/work/CT1/c1601279/rguillermin/RUN_CROCO/run_swio2_deter2_2017_2023/swio_avg_suf.nc'
grid_path = '/home/guilremy/IGE-Stochastic/Data/croco_grid_swio2.nc'
#grid_path = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/CROCO_FILES/grid/croco_grid_swio2.nc'
glorys_path = '/home/guilremy/IGE-Stochastic/Data/glorys_avg_suf.nc'
#data_path = '/lus/work/CT1/c1601279/rguillermin/GLORYS/glorys_avg_suf.nc'

lon, lat, pm, pn, msk, msk_inv, angle, _ = cplot.utils.load_grid(grid_path)
_, vorticity = cplot.utils.load_data(data_path, ('V_hor','vort'))
vorticity = vorticity[:,0,:,:]

_, vorticity_g = cplot.utils.load_data(glorys_path, ('V_hor','vort'))
vorticity_g =  vorticity_g[:,0,:,:]

lon_g, lat_g = np.meshgrid(vorticity_g.longitude, vorticity_g.latitude)

figsize = (24,12)

for index in range(vorticity.shape[0]):
    date = np.datetime_as_string(vorticity.time.data[index], 'D')
    print(f'Plotting {date}...')

    fig, axes = plt.subplots(2, 3, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    
    axes = axes.flatten()

    for ax in axes:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=11)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=11)
        ax.add_feature(cfeature.LAND, color='lightgrey', zorder=10)
        ax.set_extent(SWIO)
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), **gridline_style, zorder=10)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 8, 'color': 'k'}
        
    # =============================================================================
    # CROCO
    # =============================================================================       

        
    # =============================================================================
    # VORTICITY
    # =============================================================================
     
    levels = np.linspace(-0.00005, 0.00005, 51)
    norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)

    ax = axes[0]
    ax.set_title(f"Vorticity CROCO")
    pcm = ax.pcolormesh(lon[:, :], lat[:, :],  vorticity[index], cmap=vort_cmap, norm=norm, transform=ccrs.PlateCarree())

    cb = plt.colorbar(pcm, ax=ax, label=r'Vorticity [m²/s]', orientation='vertical')
    cb.set_ticks(cb.get_ticks())
    cb.ax.set_yticklabels(np.round(cb.get_ticks(), 6))

    # =============================================================================
    # MASKING
    # =============================================================================

    # Custom colormap: blue (negative), white (neutral), red (positive)
    mask_cmap = mcolors.ListedColormap(["blue", "white", "red"])

    # Normalization to align -1, 0, 1 to colormap indices
    mask_norm = mcolors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], mask_cmap.N)

    # Ternary mask: 0 = no mask, 1 = positive mask, -1 = negative mask
    mask_values = np.zeros_like(vorticity[index].data, dtype=int)
    mask_values[vorticity[index].data > min_vort] = 1   # Positive vorticity
    mask_values[vorticity[index].data < -min_vort] = -1 # Negative vorticity

    labeled_mask, num_features = label(mask_values)

    ax = axes[1]
    ax.set_title(rf"CROCO masked | $N_f={num_features}$")

    # Plot the mask
    pcmmask = ax.pcolormesh(lon[:, :], lat[:, :], mask_values, cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

    # Colorbar
    cb = plt.colorbar(pcmmask, ax=ax, ticks=[-1, 0, 1], orientation='vertical')
    cb.ax.set_yticklabels(['Vort < -min', '', 'Vort > min'])
    for cb_label in cb.ax.get_yticklabels():
        cb_label.set_rotation(90)
        cb_label.set_verticalalignment('center')
    cb.set_label(rf'$|$Vorticity$| > {min_vort}$')


    # =============================================================================
    # FILTERING 
    # =============================================================================

    # Initialize filtered mask as 0
    filtered_mask = np.zeros_like(mask_values, dtype=int)

    # --- Process positive vorticity regions ---
    pos_mask = (mask_values == 1)
    labeled_pos, num_pos = label(pos_mask)
    region_areas_pos = ndi_sum(pos_mask, labeled_pos, index=range(1, num_pos + 1))

    for i, area in enumerate(region_areas_pos):
        if area >= min_area:
            filtered_mask[labeled_pos == (i + 1)] = 1  # Keep as positive mask

    # --- Process negative vorticity regions ---
    neg_mask = (mask_values == -1)
    labeled_neg, num_neg = label(neg_mask)
    region_areas_neg = ndi_sum(neg_mask, labeled_neg, index=range(1, num_neg + 1))

    for i, area in enumerate(region_areas_neg):
        if area >= min_area:
            filtered_mask[labeled_neg == (i + 1)] = -1  # Keep as negative mask

    # Count total filtered features (optional)
    num_size_filtered = (np.max(labeled_pos) + np.max(labeled_neg))  # crude total count of objects before area filtering
    num_filtered_kept = len(region_areas_pos[region_areas_pos >= min_area]) + len(region_areas_neg[region_areas_neg >= min_area])

    ax = axes[2]
    ax.set_title(rf"CROCO filtered (size) | $N_f={num_filtered_kept}$")

    # Use same colormap and norm as before (centered on 0)
    pcmfilt = ax.pcolormesh(
        lon[:, :], lat[:, :], filtered_mask,
        cmap=mask_cmap, norm=mask_norm,
        transform=ccrs.PlateCarree()
    )

    # Colorbar properly labeled
    cb = plt.colorbar(pcmfilt, ax=ax, ticks=[-1, 0, 1], orientation='vertical')
    cb.ax.set_yticklabels(['Vort < -min', '', 'Vort > min'])
    for cb_label in cb.ax.get_yticklabels():
        cb_label.set_rotation(90)
        cb_label.set_verticalalignment('center')
    cb.set_label(rf'$|$Vorticity$| > {min_vort}$ & $A > {min_area}$')

    # =============================================================================
    # GLORYS
    # =============================================================================       

        
    # =============================================================================
    # VORTICITY
    # =============================================================================
     
    levels = np.linspace(-0.00005, 0.00005, 51)
    norm = mpl.colors.BoundaryNorm(levels, vort_cmap.N)

    ax = axes[3]
    ax.set_title(f"Vorticity GLORYS")
    pcm = ax.pcolormesh(lon_g[:, :], lat_g[:, :],  vorticity_g[index], cmap=vort_cmap, norm=norm, transform=ccrs.PlateCarree())

    cb = plt.colorbar(pcm, ax=ax, label=r'Vorticity [m²/s]', orientation='vertical')
    cb.set_ticks(cb.get_ticks())
    cb.ax.set_yticklabels(np.round(cb.get_ticks(), 6))

    # =============================================================================
    # MASKING
    # =============================================================================

    # Custom colormap: blue (negative), white (neutral), red (positive)
    mask_cmap = mcolors.ListedColormap(["blue", "white", "red"])

    # Normalization to align -1, 0, 1 to colormap indices
    mask_norm = mcolors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], mask_cmap.N)

    # Ternary mask: 0 = no mask, 1 = positive mask, -1 = negative mask
    mask_values = np.zeros_like(vorticity_g[index].data, dtype=int)
    mask_values[vorticity_g[index].data > min_vort] = 1   # Positive vorticity
    mask_values[vorticity_g[index].data < -min_vort] = -1 # Negative vorticity

    labeled_mask, num_features = label(mask_values)

    ax = axes[4]
    ax.set_title(rf"GLORYS masked | $N_f={num_features}$")

    # Plot the mask
    pcmmask = ax.pcolormesh(lon_g[:, :], lat_g[:, :], mask_values, cmap=mask_cmap, norm=mask_norm, transform=ccrs.PlateCarree())

    # Colorbar
    cb = plt.colorbar(pcmmask, ax=ax, ticks=[-1, 0, 1], orientation='vertical')
    cb.ax.set_yticklabels(['Vort < -min', '', 'Vort > min'])
    for cb_label in cb.ax.get_yticklabels():
        cb_label.set_rotation(90)
        cb_label.set_verticalalignment('center')
    cb.set_label(rf'$|$Vorticity$| > {min_vort}$')


    # =============================================================================
    # FILTERING 
    # =============================================================================

    # Initialize filtered mask as 0
    filtered_mask = np.zeros_like(mask_values, dtype=int)

    # --- Process positive vorticity regions ---
    pos_mask = (mask_values == 1)
    labeled_pos, num_pos = label(pos_mask)
    region_areas_pos = ndi_sum(pos_mask, labeled_pos, index=range(1, num_pos + 1))

    for i, area in enumerate(region_areas_pos):
        if area >= min_area:
            filtered_mask[labeled_pos == (i + 1)] = 1  # Keep as positive mask

    # --- Process negative vorticity regions ---
    neg_mask = (mask_values == -1)
    labeled_neg, num_neg = label(neg_mask)
    region_areas_neg = ndi_sum(neg_mask, labeled_neg, index=range(1, num_neg + 1))

    for i, area in enumerate(region_areas_neg):
        if area >= min_area:
            filtered_mask[labeled_neg == (i + 1)] = -1  # Keep as negative mask

    # Count total filtered features (optional)
    num_size_filtered = (np.max(labeled_pos) + np.max(labeled_neg))  # crude total count of objects before area filtering
    num_filtered_kept = len(region_areas_pos[region_areas_pos >= min_area]) + len(region_areas_neg[region_areas_neg >= min_area])

    ax = axes[5]
    ax.set_title(rf"GLORYS filtered (size) | $N_f={num_filtered_kept}$")

    # Use same colormap and norm as before (centered on 0)
    pcmfilt = ax.pcolormesh(
        lon_g[:, :], lat_g[:, :], filtered_mask,
        cmap=mask_cmap, norm=mask_norm,
        transform=ccrs.PlateCarree()
    )

    # Colorbar properly labeled
    cb = plt.colorbar(pcmfilt, ax=ax, ticks=[-1, 0, 1], orientation='vertical')
    cb.ax.set_yticklabels(['Vort < -min', '', 'Vort > min'])
    for cb_label in cb.ax.get_yticklabels():
        cb_label.set_rotation(90)
        cb_label.set_verticalalignment('center')
    cb.set_label(rf'$|$Vorticity$| > {min_vort}$ & $A > {min_area}$')

    fig.suptitle(f'SWIO {date}')
    plt.tight_layout()
    
    cplot.utils.save_figure(fig, f"frame_{date}", isFilm=True, filmDir='vorticity_filtered')
    plt.close()
