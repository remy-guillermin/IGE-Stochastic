import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from .utils import load_grid, load_data, save_figure

def eke(data_files, boxes, names, colors, y_min=1e2, y_max=1e4, vertical_y_min=1e-3, vertical_y_max=1e-1):
    """
    Calculate and plot the time series of EKE for specified regions.

    Parameters
    ----------
    data_path : str
        Path to the simulation data directory.
    simu : str
        Simulation subdirectory.
    grid_path : str
        Path to the grid file.
    boxes : list of tuplesda
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    """
    # Load grid data
    lon, lat, _, _, msk, msk_inv, angle, h = load_grid()

    eke_results = {name: [] for name in names}
    vertical_eke_results = {name: [] for name in names}
    time_results = []

    for data_file in data_files:
        # Load simulation data
        u, v, w, time, s_rho = load_data(data_file, ('u', 'v', 'w', 'time', 's_rho'))

        depth = h * s_rho # Profondeur
        depth = np.transpose(depth.data, (2, 0, 1)) # Transpose depth to match u, v, w shape
        ddepth = np.diff(depth, axis=0)
        
        fill_value = 9.96921e+36
        u = u.where((u != fill_value), np.nan)
        v = v.where((v != fill_value), np.nan)
        w = w.where((w != fill_value), np.nan)

        # Moyenne annuelle
        u_yr = np.mean(u, axis = 0)
        v_yr = np.mean(v, axis = 0)
        w_yr = np.mean(w, axis = 0)

        # Vitesse turbulente
        ut = u_yr - u
        vt = v_yr - v
        wt = w_yr - w

        # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
        angle_expand = angle[:,:].data.reshape(1, angle.shape[0], angle.shape[1], 1)
        
        ut_geo = ut[:,:-1,:,:].data * np.cos(angle_expand[:,:-1,:-1:,]) - vt[:,:,:-1,:].data * np.sin(angle_expand[:,:-1,:-1,:])
        vt_geo = ut[:,:-1,:,:].data * np.sin(angle_expand[:,:-1,:-1:,]) + vt[:,:,:-1,:].data * np.cos(angle_expand[:,:-1,:-1,:])
        wt_geo = wt[:,:-1,:-1,:].data

        # Calculating EKE
        EKE = 1 / 2 * (ut_geo ** 2 + vt_geo ** 2 + wt_geo ** 2)
        
        # Integrating EKE over depth
        ddepth_expand = ddepth[:,:-1,:-1].reshape(EKE.shape[0]-1, EKE.shape[1], EKE.shape[2], 1)
        EKE_weighted = EKE[1:,:,:] * ddepth_expand
        EKE_integrated = np.sum(EKE_weighted, axis=0)

        # Extracting EKE for each box
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))[:-1,:-1]
            EKE_box = EKE[:,box_mask,:]
            ddepth_box = ddepth[:,:-1,:-1][:,box_mask]
            EKE_integrated_box = EKE_integrated[box_mask]

            # Calculate the spatial mean of EKE_box over time
            EKE_box_sum = np.nansum(EKE_box, axis=1)
            EKE_integrated_box_sum = np.nansum(EKE_integrated_box, axis=0)
            eke_results[name].append(EKE_box_sum.T)
            vertical_eke_results[name].append(EKE_integrated_box_sum.T)

        time_results.append(time.data)

    # Concatenate results
    for name in names:
        eke_results[name] = np.concatenate(eke_results[name])
        vertical_eke_results[name] = np.concatenate(vertical_eke_results[name])
    time_results = np.concatenate(time_results)

    # Plot EKE time series
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.semilogy(time_results, eke_results[name][:,-1], label=name, color=color)
        ax.set_ylabel('EKE [m²/s²]')
        ax.legend()
        ax.set_ylim(y_min, y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('EKE Over Time for Different Boxes')
    save_figure(fig, f"eke_time_series.png")
    plt.close(fig)

    # Plot vertical EKE time series
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.semilogy(time_results, vertical_eke_results[name], label=name, color=color)
        ax.set_ylabel('Vertically integrated EKE [m²/s²]')
        ax.legend()
        ax.set_ylim(vertical_y_min, vertical_y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Vertically integrated EKE Over Time for Different Boxes')
    save_figure(fig, f"vertical_eke_time_series.png")
    plt.close(fig)
