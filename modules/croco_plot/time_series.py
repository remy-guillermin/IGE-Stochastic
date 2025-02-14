"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage des données temporelles de CROCO.
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from .utils import load_grid, load_data, save_figure

def box_eke(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        y_min=1e2, 
        y_max=1e4, 
        vertical_y_min=1e-3, 
        vertical_y_max=1e-1):
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
    if isinstance(data_files, str):
        data_files = [data_files]
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
        ax.set_ylabel('EKE [$m^2/s^2$]')
        ax.legend()
        ax.set_ylim(y_min, y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('EKE Over Time for Different Boxes')
    save_figure(fig, f"eke_boxes_time_series.png")
    plt.close(fig)

    # Plot vertical EKE time series
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.semilogy(time_results, vertical_eke_results[name], label=name, color=color)
        ax.set_ylabel('[$m^3/s^2$]')
        ax.legend()
        ax.set_ylim(vertical_y_min, vertical_y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Vertically integrated EKE Over Time for Different Boxes')
    save_figure(fig, f"vertical_eke_boxes_time_series.png")
    plt.close(fig)
    
def box_mke(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        y_min=1e2, 
        y_max=1e4, 
        vertical_y_min=1e5, 
        vertical_y_max=1e7):
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
    if isinstance(data_files, str):
        data_files = [data_files]
    # Load grid data
    lon, lat, _, _, msk, msk_inv, angle, h = load_grid()

    mke_results = {name: [] for name in names}
    vertical_mke_results = {name: [] for name in names}
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

        # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
        angle_expand = angle[:,:].data.reshape(1, 1, angle.shape[0], angle.shape[1])
        
        u_geo = u[:,:,:-1,:].data * np.cos(angle_expand[:,:,:-1:,:-1]) - v[:,:,:,:-1].data * np.sin(angle_expand[:,:,:-1,:-1])
        v_geo = u[:,:,:-1,:].data * np.sin(angle_expand[:,:,:-1:,:-1]) + v[:,:,:,:-1].data * np.cos(angle_expand[:,:,:-1,:-1])
        w_geo = w[:,:,:-1,:-1].data

        # Calculating EKE
        MKE = 1 / 2 * (u_geo ** 2 + v_geo ** 2 + w_geo ** 2)

        # Extracting EKE for each box
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))[:-1,:-1]
            MKE_box = MKE[:,:,box_mask]
            ddepth_box = ddepth[:,:-1,:-1][:,box_mask]
            MKE_weighted = MKE_box[:,1:,:] * ddepth_box
            MKE_integrated = np.sum(MKE_weighted, axis=1)

            # Calculate the spatial mean of EKE_box over time
            MKE_box_sum = np.nansum(MKE_box[:,-1,:], axis=1)
            MKE_integrated_box_sum = np.nansum(MKE_integrated, axis=1)
            mke_results[name].append(MKE_box_sum)
            vertical_mke_results[name].append(MKE_integrated_box_sum)

        time_results.append(time.data)

    # Concatenate results
    for name in names:
        mke_results[name] = np.concatenate(mke_results[name])
        vertical_mke_results[name] = np.concatenate(vertical_mke_results[name])
    time_results = np.concatenate(time_results)

    # Plot EKE time series
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.semilogy(time_results, mke_results[name][:], label=name, color=color)
        ax.set_ylabel('MKE [$m^2/s^2$]')
        ax.legend()
        ax.set_ylim(y_min, y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('MKE Over Time for Different Boxes')
    save_figure(fig, f"mke_boxes_time_series.png")
    plt.close(fig)

    # Plot vertical EKE time series
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.semilogy(time_results, vertical_mke_results[name], label=name, color=color)
        ax.set_ylabel('[$m^3/s^2$]')
        ax.legend()
        ax.set_ylim(vertical_y_min, vertical_y_max)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Vertically integrated EKE Over Time for Different Boxes')
    save_figure(fig, f"vertical_mke_boxes_time_series.png")
    plt.close(fig)

def box_sla(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Level Anomaly (SLA) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    sla_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        zeta, time = load_data(data_file, ('zeta', 'time'))
        
        fill_value = 9.96921e+36
        zeta = zeta.where((zeta != fill_value), np.nan)
        
        zeta_yr = np.nanmean(zeta, axis=0)
        
        SLA = np.array(zeta - zeta_yr)
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            SLA_box = SLA[:,box_mask]

            # Calculate the spatial mean of SLA_box over time
            SLA_box_mean = np.nanmean(SLA_box, axis=1)
            sla_results[name].append(SLA_box_mean)

        time_results.append(time.data)
        
    for name in names:
        sla_results[name] = np.concatenate(sla_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, sla_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean
        sla_rolling_mean = np.convolve(sla_results[name], np.ones(roll)/roll, mode='same')
        valid_indices = ~np.isnan(sla_rolling_mean)
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sla_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('SLA [m]')
        ax.set_ylim(y_min, y_max)
        ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Level Anomaly Over Time for Different Boxes')
    save_figure(fig, f"sla_boxes_time_series.png")
    plt.close(fig)
        
def box_sta(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Temperature Anomaly (STA) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    sta_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        temp, time = load_data(data_file, ('temp', 'time'))
        temp = temp[:,-1,:,:]
        
        fill_value = 9.96921e+36
        temp = temp.where((temp != fill_value), np.nan)
        
        temp_yr = np.nanmean(temp, axis=0)
        
        STA = np.array(temp - temp_yr)
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            STA_box = STA[:,box_mask]

            # Calculate the spatial mean of STA_box over time
            STA_box_mean = np.nanmean(STA_box, axis=1)
            sta_results[name].append(STA_box_mean)

        time_results.append(time.data)
        
    for name in names:
        sta_results[name] = np.concatenate(sta_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, sta_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean
        sta_rolling_mean = np.convolve(sta_results[name], np.ones(roll)/roll, mode='same')
        valid_indices = ~np.isnan(sta_rolling_mean)
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sta_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('STA [°C]')
        ax.set_ylim(y_min, y_max)
        ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Temperature Anomaly Over Time for Different Boxes')
    save_figure(fig, f"sta_boxes_time_series.png")
    plt.close(fig)
    
def box_ssa(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Salinity Anomaly (SSA) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    ssa_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        salt, time = load_data(data_file, ('salt', 'time'))
        salt = salt[:,-1,:,:]
        
        fill_value = 9.96921e+36
        salt = salt.where((salt != fill_value), np.nan)
        
        salt_yr = np.nanmean(salt, axis=0)
        
        SSA = np.array(salt - salt_yr)
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            SSA_box = SSA[:,box_mask]

            # Calculate the spatial mean of SSA_box over time
            SSA_box_mean = np.nanmean(SSA_box, axis=1)
            ssa_results[name].append(SSA_box_mean)

        time_results.append(time.data)
        
    for name in names:
        ssa_results[name] = np.concatenate(ssa_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, ssa_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean
        ssa_rolling_mean = np.convolve(ssa_results[name], np.ones(roll)/roll, mode='same')
        valid_indices = ~np.isnan(ssa_rolling_mean)
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], ssa_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('SSA [psu]')
        ax.set_ylim(y_min, y_max)
        ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Salinity Anomaly Over Time for Different Boxes')
    save_figure(fig, f"ssa_boxes_time_series.png")
    plt.close(fig)
    
def box_ssh(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Surface Height (SSH) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    ssh_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        zeta, time = load_data(data_file, ('zeta', 'time'))
        
        fill_value = 9.96921e+36
        zeta = zeta.where((zeta != fill_value), np.nan).data
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            ssh_box = zeta[:,box_mask]

            # Calculate the spatial mean of SSH over time
            ssh_box_mean = np.nanmean(ssh_box, axis=1)
            ssh_results[name].append(ssh_box_mean)

        time_results.append(time.data)
        
    for name in names:
        ssh_results[name] = np.concatenate(ssh_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, ssh_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean
        ssh_rolling_mean = np.convolve(ssh_results[name], np.ones(roll)/roll, mode='same')
        valid_indices = ~np.isnan(ssh_rolling_mean)
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], ssh_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('SSH [m]')
        ax.set_ylim(y_min, y_max)
        ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Surface Height Over Time for Different Boxes')
    save_figure(fig, f"ssh_boxes_time_series.png")
    plt.close(fig)
        
def box_sst(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Surface Temperature (SST) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    sst_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        temp, time = load_data(data_file, ('temp', 'time'))
        temp = temp[:,-1,:,:]
        
        fill_value = 9.96921e+36
        temp = temp.where((temp != fill_value), np.nan).data
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            sst_box = temp[:,box_mask]

            # Calculate the spatial mean of STA_box over time
            sst_box_mean = np.nanmean(sst_box, axis=1)
            sst_results[name].append(sst_box_mean)

        time_results.append(time.data)
        
    for name in names:
        sst_results[name] = np.concatenate(sst_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, sst_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean 
        sst_rolling_mean = np.convolve(sst_results[name], np.ones(roll)/roll, mode='same')
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sst_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('SST [°C]')
        ax.set_ylim(y_min, y_max)
        #ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Surface Temperature Over Time for Different Boxes')
    save_figure(fig, f"sta_boxes_time_series.png")
    plt.close(fig)
    
def box_sss(data_files, 
        boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
        names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
        colors=['saddlebrown', 'darkorchid', 'navy', 'teal'], 
        roll = 9,
        y_min=None, 
        y_max=None):
    """
    Calculate and plot the time series of Sea Surface Salinity (SSS) for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    boxes : list of tuples
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list of str
        Names of the regions.
    colors : list of str
        Colors for the plot lines.
    roll : int
        Window size for the centered rolling mean.
    y_min : float
        Minimum y-axis value for the plot.
    y_max : float
        Maximum y-axis value for the plot.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
    lon, lat, _, _, msk, msk_inv, _, _ = load_grid()
    
    sss_results = {name: [] for name in names}
    time_results = []
    
    for data_file in data_files:
        # Load simulation data
        salt, time = load_data(data_file, ('salt', 'time'))
        salt = salt[:,-1,:,:]
        
        fill_value = 9.96921e+36
        salt = salt.where((salt != fill_value), np.nan).data
        
        for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            sss_box = salt[:,box_mask]

            # Calculate the spatial mean of SSA_box over time
            sss_box_mean = np.nanmean(sss_box, axis=1)
            sss_results[name].append(sss_box_mean)

        time_results.append(time.data)
        
    for name in names:
        sss_results[name] = np.concatenate(sss_results[name])
    time_results = np.concatenate(time_results)
    
    fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

    for ax, (name, color) in zip(axes, zip(names, colors)):
        ax.plot(time_results, sss_results[name], color=color, linestyle='--', linewidth=1)
        # Apply centered rolling mean
        sss_rolling_mean = np.convolve(sss_results[name], np.ones(roll)/roll, mode='same')
        ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sss_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
        ax.set_ylabel('SSS [psu]')
        ax.set_ylim(y_min, y_max)
        ax.set_title(name)

    axes[-1].set_xlabel('Time')
    fig.suptitle('Sea Surface Salinity Over Time for Different Boxes')
    save_figure(fig, f"sss_boxes_time_series.png")
    plt.close(fig)
