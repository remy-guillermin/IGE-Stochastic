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

def multiple_time_series(data_files, 
                         variables='all',
                         boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
                         names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
                         colors=['saddlebrown', 'darkorchid', 'navy', 'teal'],
                         roll = 9):
    """
    Calculate and plot multiple time series for specified regions.

    Parameters
    ----------
    data_files : list of str
        List of paths to the simulation data files.
    variables : list of str, optional
        List of variables to plot, by default 'all'
        options : 'ke', 'eke', 'mke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst'
    boxes : list, optional
        List of tuples defining the regions (lon1, lon2, lat1, lat2).
    names : list, optional
       Names of the regions.
    colors : list, optional
        Colors for the plot lines.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
        
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['ke', 'eke', 'mke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']
        
    time_results = []
    if 'ke' in variables:
        ke_results = []
        vertical_ke_results = []
    
    if 'eke' in variables:
        eke_results = {name: [] for name in names}
        vertical_eke_results = {name: [] for name in names}
    
    if 'mke' in variables:
        mke_results = {name: [] for name in names}
        vertical_mke_results = {name: [] for name in names}
    
    if 'sla' in variables:
        sla_results = {name: [] for name in names}
    
    if 'ssa' in variables:
        ssa_results = {name: [] for name in names}
    
    if 'sta' in variables:
        sta_results = {name: [] for name in names}
    
    if 'ssh' in variables:
        ssh_results = {name: [] for name in names}
    
    if 'sss' in variables:
        sss_results = {name: [] for name in names}
    
    if 'sst' in variables:
        sst_results = {name: [] for name in names}
        
    # Load grid data
    lon, lat, _, _, _, _, _, h = load_grid()
    
    fill_value = 9.96921e+36
    
    for data_file in data_files:
        zeta, temp, salt, u, v, w, time, s_rho = load_data(data_file, ('zeta', 'temp', 'salt', 'u', 'v', 'w', 'time', 's_rho'))
        if 'ke' in variables or 'eke' in variables or 'mke' in variables:
            print('Converting velocities')
            u, v, w = u.data, v.data, w.data
            u[u == fill_value] = np.nan
            v[v == fill_value] = np.nan
            w[w == fill_value] = np.nan
        if 'sla' in variables or 'ssh' in variables:
            print('Converting free surface')
            zeta = zeta.data
            zeta[zeta == fill_value] = np.nan
        if 'sta' in variables or 'sst' in variables:
            print('Converting temperature')
            temp = temp.data
            temp[temp == fill_value] = np.nan
            temp = temp[:,-1,:,:]
        if 'ssa' in variables or 'sss' in variables:
            print('Converting salinity')
            salt = salt.data
            salt[salt == fill_value] = np.nan
            salt = salt[:,-1,:,:]
        
        time_results.append(time.data)
        time = None
        
        depth = h * s_rho # Profondeur
        depth = np.transpose(depth.data, (2, 0, 1)) # Transpose depth to match u, v, w shape
        ddepth = np.diff(depth, axis=0)
        ddepth = ddepth.reshape(1,ddepth.shape[0], ddepth.shape[1], ddepth.shape[2])
        depth, s_rho = None, None
        
        if 'ke' in variables:
            print('Calculating KE')
            KE = 1 / 2 * (u[:,:,:-1,:] ** 2 + v[:,:,:,:-1] ** 2 + w[:,:,:-1,:-1] ** 2)
            surf_KE = np.nansum(KE[:,-1,:,:], axis=(1, 2))
            KE_weighted = KE[:,1:,:,:] * ddepth[:,:,:-1,:-1]
            KE_integrated = np.nansum(KE_weighted, axis=1)
            ke_results.append(surf_KE)
            vertical_ke_results.append(np.nansum(KE_integrated, axis=(1,2)))
            KE, surf_KE, KE_weighted, KE_integrated = None, None, None, None
            
        if 'mke' in variables:
            print('Calculating MKE')
            MKE = 1 / 2 * (u[:,:,:-1,:] ** 2 + v[:,:,:,:-1] ** 2 + w[:,:,:-1,:-1] ** 2)
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))[:-1,:-1]
                MKE_box = MKE[:,:,box_mask]
                ddepth_box = ddepth[:,:,:-1,:-1][:,:,box_mask]
                MKE_weighted = MKE_box[:,1:,:] * ddepth_box
                MKE_integrated = np.sum(MKE_weighted, axis=1)
                MKE_box_sum = np.nansum(MKE_box[:,-1,:], axis=1)
                MKE_integrated_box_sum = np.nansum(MKE_integrated, axis=1)
                mke_results[name].append(MKE_box_sum)
                vertical_mke_results[name].append(MKE_integrated_box_sum)
            MKE, MKE_box, MKE_weighted, MKE_integrated = None, None, None, None
            
        if 'eke' in variables:
            print('Calculating EKE')
            ut = (u - np.nanmean(u, axis=0))
            vt = (v - np.nanmean(v, axis=0))
            wt = (w - np.nanmean(w, axis=0))
            EKE = 1 / 2 * (ut[:,:,:-1,:] ** 2 + vt[:,:,:,:-1] ** 2 + wt[:,:,:-1,:-1] ** 2)
            ut, vt, wt = None, None, None
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))[:-1,:-1]
                EKE_box = EKE[:,:,box_mask]
                ddepth_box = ddepth[:,:,:-1,:-1][:,:,box_mask]
                EKE_weighted = EKE_box[:,1:,:] * ddepth_box
                EKE_integrated = np.nansum(EKE_weighted, axis=1)
                EKE_box_sum = np.nansum(EKE_box[:,-1,:], axis=1)
                EKE_integrated_box_sum = np.nansum(EKE_integrated, axis=1)
                eke_results[name].append(EKE_box_sum)
                vertical_eke_results[name].append(EKE_integrated_box_sum)
            EKE, EKE_box, EKE_weighted, EKE_integrated = None, None, None, None
            
        if 'sla' in variables:
            print('Calculating SLA')
            SLA = (zeta - np.nanmean(zeta, axis=0))
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                SLA_box = SLA[:,box_mask]
                SLA_box_mean = np.nanmean(SLA_box, axis=1)
                sla_results[name].append(SLA_box_mean)
            SLA = None
            
        if 'ssa' in variables:
            print('Calculating SSA')
            SSA = (salt - np.nanmean(salt, axis=0))
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                SSA_box = SSA[:,box_mask]
                SSA_box_mean = np.nanmean(SSA_box, axis=1)
                ssa_results[name].append(SSA_box_mean)
            SSA = None
            
        if 'sta' in variables:
            print('Calculating STA')
            STA = (temp - np.nanmean(temp, axis=0))
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                STA_box = STA[:,box_mask]
                STA_box_mean = np.nanmean(STA_box, axis=1)
                sta_results[name].append(STA_box_mean)
            STA = None
            
        if 'ssh' in variables:
            print('Calculating SSH')
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                SSH_box = zeta[:,box_mask]
                SSH_box_mean = np.nanmean(SSH_box, axis=1)
                ssh_results[name].append(SSH_box_mean)
        
        if 'sss' in variables:
            print('Calculating SSS')
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                SSS_box = salt[:,box_mask]
                SSS_box_mean = np.nanmean(SSS_box, axis=1)
                sss_results[name].append(SSS_box_mean)
        
        if 'sst' in variables:
            print('Calculating SST')
            for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                SST_box = temp[:,box_mask]
                SST_box_mean = np.nanmean(SST_box, axis=1)
                sst_results[name].append(SST_box_mean)
    
    time_results = np.concatenate(time_results)
    
    if 'ke' in variables:
        ke_results = np.concatenate(ke_results)
        vertical_ke_results = np.concatenate(vertical_ke_results)
        
        fig, axes = plt.subplots(2, 1, figsize = (12, 4), sharex= True)
        
        ax = axes[0]
        ax.semilogy(time_results, ke_results, color='black')
        ax.set_ylabel('KE [$m^2/s^2$]')
        
        ax = axes[1]
        ax.semilogy(time_results, vertical_ke_results, color='black')
        ax.set_ylabel('[$m^3/s^2$]')
        
        axes[-1].set_xlabel('Time')
        fig.suptitle('KE Over Time for the Whole Domain')
        save_figure(fig, f"ke_global_time_series.png")
        plt.close(fig)
        
    if 'mke' in variables:
        for name in names:
            mke_results[name] = np.concatenate(mke_results[name])
            vertical_mke_results[name] = np.concatenate(vertical_mke_results[name])
        
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)
        
        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.semilogy(time_results, mke_results[name][:], label=name, color=color)
            ax.set_ylabel('MKE [$m^2/s^2$]')
            ax.legend()
        
        axes[-1].set_xlabel('Time')
        fig.suptitle('MKE Over Time for Different Boxes')
        save_figure(fig, f"mke_boxes_time_series.png")
        plt.close(fig)
        
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)
        
        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.semilogy(time_results, vertical_mke_results[name], label=name, color=color)
            ax.set_ylabel('[$m^3/s^2$]')
            ax.legend()
        
        axes[-1].set_xlabel('Time')
        fig.suptitle('Vertically integrated MKE Over Time for Different Boxes')
        save_figure(fig, f"vertical_mke_boxes_time_series.png")
        plt.close(fig)
    
    if 'eke' in variables:
        # Concatenate results
        for name in names:
            eke_results[name] = np.concatenate(eke_results[name])
            vertical_eke_results[name] = np.concatenate(vertical_eke_results[name])

        # Plot EKE time series
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.semilogy(time_results, eke_results[name][:], label=name, color=color)
            ax.set_ylabel('MKE [$m^2/s^2$]')
            ax.legend()

        axes[-1].set_xlabel('Time')
        fig.suptitle('MKE Over Time for Different Boxes')
        save_figure(fig, f"mke_boxes_time_series.png")
        plt.close(fig)

        # Plot vertical EKE time series
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.semilogy(time_results, vertical_eke_results[name], label=name, color=color)
            ax.set_ylabel('[$m^3/s^2$]')
            ax.legend()

        axes[-1].set_xlabel('Time')
        fig.suptitle('Vertically integrated EKE Over Time for Different Boxes')
        save_figure(fig, f"vertical_mke_boxes_time_series.png")
        plt.close(fig)
        
    if 'sla' in variables:
        for name in names:
            sla_results[name] = np.concatenate(sla_results[name])
        
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)
        
        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, sla_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean
            sla_rolling_mean = np.convolve(sla_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sla_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('SLA [m]')
            ax.set_title(name)
        
        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Level Anomaly Over Time for Different Boxes')
        save_figure(fig, f"sla_boxes_time_series.png")
        plt.close(fig)
        
    if 'sta' in variables:
        for name in names:
            sta_results[name] = np.concatenate(sta_results[name])
        
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, sta_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean
            sta_rolling_mean = np.convolve(sta_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sta_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('STA [°C]')
            ax.set_title(name)

        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Temperature Anomaly Over Time for Different Boxes')
        save_figure(fig, f"sta_boxes_time_series.png")
        plt.close(fig)
        
    if 'ssa' in variables:
        for name in names:
            ssa_results[name] = np.concatenate(ssa_results[name])

        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, ssa_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean
            ssa_rolling_mean = np.convolve(ssa_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], ssa_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('SSA [psu]')
            ax.set_title(name)

        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Salinity Anomaly Over Time for Different Boxes')
        save_figure(fig, f"ssa_boxes_time_series.png")
        plt.close(fig)
        
    if 'ssh' in variables:
        for name in names:
            ssh_results[name] = np.concatenate(ssh_results[name])

        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, ssh_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean
            ssh_rolling_mean = np.convolve(ssh_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], ssh_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('SSH [m]')
            ax.set_title(name)

        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Surface Height Over Time for Different Boxes')
        save_figure(fig, f"ssh_boxes_time_series.png")
        plt.close(fig)
    
    if 'sss' in variables:
        for name in names:
            sss_results[name] = np.concatenate(sss_results[name])

        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, sss_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean
            sss_rolling_mean = np.convolve(sss_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sss_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('SSS [psu]')
            ax.set_title(name)

        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Surface Salinity Over Time for Different Boxes')
        save_figure(fig, f"sss_boxes_time_series.png")
        plt.close(fig)
    
    if 'sst' in variables:
        for name in names:
            sst_results[name] = np.concatenate(sst_results[name])
        
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

        for ax, (name, color) in zip(axes, zip(names, colors)):
            ax.plot(time_results, sst_results[name], color=color, linestyle='--', linewidth=1)
            # Apply centered rolling mean 
            sst_rolling_mean = np.convolve(sst_results[name], np.ones(roll)/roll, mode='same')
            ax.plot(time_results[int((roll-1)/2):-int((roll-1)/2)], sst_rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=1.5, alpha=0.8)
            ax.set_ylabel('SST [°C]')
            ax.set_title(name)

        axes[-1].set_xlabel('Time')
        fig.suptitle('Sea Surface Temperature Over Time for Different Boxes')
        save_figure(fig, f"sst_boxes_time_series.png")
        plt.close(fig)