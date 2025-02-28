"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage des données temporelles de CROCO.
"""

import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from .utils import load_grid, load_data, save_figure, plot_time_series

def multiple_time_series(data_files, 
                         variables='all',
                         boxes=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
                         names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
                         colors=['saddlebrown', 'darkorchid', 'navy', 'teal'],
                         roll = 9,
                         grid_path=None,
                         interactive=False):
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
    grid_path : str, optional
        Path to the grid file, by default None
    interactive : bool, optional
        Whether to use an interactive backend for plotting.
    """
    if isinstance(data_files, str):
        data_files = [data_files]
        
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['ke', 'eke', 'mke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']
        
    time_results = []
    results = {var: {name: [] for name in names} for var in variables if var != 'ke'}
    ke_results = []

    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, _, _, h = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, _, _, h = load_grid()
    
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
        
        cell_surface = 1 / (pm * pn).data
        cell_surface[(1-msk).astype(int)] = np.nan
        
        depth = h * s_rho # Profondeur
        depth = np.transpose(depth.data, (2, 0, 1)) # Transpose depth to match u, v, w shape
        ddepth = np.diff(depth, axis=0)
        ddepth = ddepth.reshape(1,ddepth.shape[0], ddepth.shape[1], ddepth.shape[2])
        cell_volume = ddepth * cell_surface
        domain_volume = np.nansum(cell_volume)
        depth, ddepth, cell_surface, s_rho = None, None, None, None
        
        if 'ke' in variables:
            print('Calculating KE')
            KE = 1 / 2 * (u[:,:,:-1,:] ** 2 + v[:,:,:,:-1] ** 2 + w[:,:,:-1,:-1] ** 2)
            surf_KE = np.nansum(KE[:,-1,:,:], axis=(1, 2))
            KE_weighted = KE[:,1:,:,:] * cell_volume[:,:,:-1,:-1] / domain_volume
            ke_results.append(np.sum(KE_weighted, axis=(1, 2, 3)))
            del KE, surf_KE, KE_weighted
        
        for var in ['mke', 'eke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']:
            if var in variables:
                print(f'Calculating {var.upper()}')
                if var == 'mke':
                    var_data = 1 / 2 * (u[:,:,:-1,:] ** 2 + v[:,:,:,:-1] ** 2 + w[:,:,:-1,:-1] ** 2)
                elif var == 'eke':
                    ut = (u - np.nanmean(u, axis=0))
                    vt = (v - np.nanmean(v, axis=0))
                    wt = (w - np.nanmean(w, axis=0))
                    var_data = 1 / 2 * (ut[:,:,:-1,:] ** 2 + vt[:,:,:,:-1] ** 2 + wt[:,:,:-1,:-1] ** 2)
                    del ut, vt, wt
                elif var == 'sla':
                    var_data = (zeta - np.nanmean(zeta, axis=0))
                elif var == 'ssa':
                    var_data = (salt - np.nanmean(salt, axis=0))
                elif var == 'sta':
                    var_data = (temp - np.nanmean(temp, axis=0))
                elif var == 'ssh':
                    var_data = zeta
                elif var == 'sss':
                    var_data = salt
                elif var == 'sst':
                    var_data = temp
                
                for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
                    box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                    if var in ['mke', 'eke']:
                        box_mask = box_mask[:-1,:-1]
                    if var in ['sla', 'ssh', 'sta', 'ssa', 'sss', 'sst']:
                        box_data = var_data[:,box_mask]
                    else:
                        box_data = var_data[:,:,box_mask]
                    if var in ['mke', 'eke']:
                        volume_box = cell_volume[:,:,:-1,:-1][:,:,box_mask]
                        weighted_data = box_data[:,1:,:] * volume_box / domain_volume
                        box_sum = np.nansum(weighted_data[:,-1,:], axis=1)
                    else:
                        box_sum = np.nanmean(box_data, axis=1)
                        weighted_data = None
                    results[var][name].append(box_sum)
                del var_data, box_data, weighted_data

    time_results = np.concatenate(time_results)
    
    if 'ke' in variables:
        ke_results = np.concatenate(ke_results)
        fig, ax = plt.subplots(1, 1, figsize=(12, 4))
        ax.semilogy(time_results, ke_results, color='black')
        ax.set_ylabel('KE [$m^2/s^2$]')
        ax.set_xlabel('Time')
        fig.suptitle('KE Over Time for the Whole Domain')
        save_figure(fig, f"ke_global_time_series.png")
        plt.close(fig)
    
    for var, ylabel, filename, title in [
        ('mke', 'MKE [$m^2/s^2$]', 'mke_boxes_time_series.png', 'MKE Over Time for Different Boxes'),
        ('eke', 'EKE [$m^2/s^2$]', 'eke_boxes_time_series.png', 'EKE Over Time for Different Boxes'),
        ('sla', 'SLA [m]', 'sla_boxes_time_series.png', 'Sea Level Anomaly Over Time for Different Boxes'),
        ('sta', 'STA [°C]', 'sta_boxes_time_series.png', 'Sea Temperature Anomaly Over Time for Different Boxes'),
        ('ssa', 'SSA [psu]', 'ssa_boxes_time_series.png', 'Sea Salinity Anomaly Over Time for Different Boxes'),
        ('ssh', 'SSH [m]', 'ssh_boxes_time_series.png', 'Sea Surface Height Over Time for Different Boxes'),
        ('sss', 'SSS [psu]', 'sss_boxes_time_series.png', 'Sea Surface Salinity Over Time for Different Boxes'),
        ('sst', 'SST [°C]', 'sst_boxes_time_series.png', 'Sea Surface Temperature Over Time for Different Boxes')
    ]:
        if var in variables:
            for name in names:
                results[var][name] = np.concatenate(results[var][name])
            plot_time_series(time_results, results[var], ylabel, roll, names, colors, filename, title, interactive=interactive)
            
def time_series_from_dataset(
    datasets,
    variables='all',
    figwidth=12,
    names=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'],
    roll = 9,
    interactive: bool = False
):   
    if interactive:
        original_backend = matplotlib.get_backend()  # Store the original backend
        matplotlib.use('tkagg')  # Temporarily switch to interactive backend
           
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['Energies', 'Anomalies', 'Fields']
    
    for name in names:
        files = datasets[name]
        if 'Energies' in variables:
            fig_energy, axes_energy = plt.subplots(2, 1, figsize=(figwidth, int(2*figwidth/3)), sharex=True)
        if 'Anomalies' in variables:
            fig_anomaly, axes_anomaly = plt.subplots(3, 1, figsize=(figwidth, figwidth), sharex=True)
        if 'Fields' in variables:
            fig_field, axes_field = plt.subplots(3, 1, figsize=(figwidth, figwidth), sharex=True)

        for file in files:
            source = Path(file).parent.name
            print(source)
            f = xr.open_dataset(file)
            time = pd.to_datetime(f.time.data)
            
            if 'mercator' in file:
                style = 'r.'
            else:
                style = '-'
            if 'Energies' in variables:
                ax = axes_energy[0]
                mke = f.mke.data
                
                ax.semilogy(time, mke, style, label=source)
                ax.set_title('Mean Kinetic Energy')
                ax.set_ylabel('MKE [$m^2/s^2$]')
                ax.set_xlabel('Time')
                
                ax = axes_energy[1]
                eke = f.eke.data
                
                ax.semilogy(time, eke, style, label=source)
                ax.set_title('Eddy Kinetic Energy')
                ax.set_ylabel('EKE [$m^2/s^2$]')
                ax.set_xlabel('Time')
            if 'Anomalies' in variables:
                ax = axes_anomaly[0]
                sla = f.sla.data
                ax.plot(time, sla, style, label=source)
                ax.set_title('Sea Level Anomaly')
                ax.set_ylabel('SLA [$m$]')
                ax.set_xlabel('Time')
                
                ax = axes_anomaly[1]
                ssa = f.ssa.data
                ax.plot(time, ssa, style, label=source)
                ax.set_title('Sea Salinity Anomaly')
                ax.set_ylabel('SSA [$psu$]')
                ax.set_xlabel('Time')
                
                ax = axes_anomaly[2]
                sta = f.sta.data
                ax.plot(time, sta, style, label=source)
                ax.set_title('Sea Temperature Anomaly')
                ax.set_ylabel('STA [$°C$]')
                ax.set_xlabel('Time')
            if 'Fields' in variables:
                ax = axes_field[0]
                ssh = f.ssh.data
                
                ax.plot(time, ssh, style, label=source)
                ax.set_title('Sea Surface Height')
                ax.set_ylabel('SSH [$m$]')
                ax.set_xlabel('Time')
                
                ax = axes_field[1]
                sss = f.sss.data
                
                ax.plot(time, sss, style, label=source)
                ax.set_title('Sea Surface Salinity')
                ax.set_ylabel('SSS [$psu$]')
                ax.set_xlabel('Time')
                
                ax = axes_field[2]
                sst = f.sst.data
                
                ax.plot(time, sst, style, label=source)
                ax.set_title('Sea Surface Temperature')
                ax.set_ylabel('SST [$°C$]')
                ax.set_xlabel('Time')    
            f.close()
        
        fig_energy.suptitle(f'Energies over time for {name}')
        fig_anomaly.suptitle(f'Anomalies over time for {name}')
        fig_field.suptitle(f'Fields over time for {name}')
        
        axes_energy[0].legend(loc="upper right")
        axes_anomaly[0].legend(loc="upper right")
        axes_field[0].legend(loc="upper right")
        
        fig_energy.tight_layout()
        fig_anomaly.tight_layout()
        fig_field.tight_layout()
        
        save_figure(fig_energy, f"{name}_energies_time_series.png")
        save_figure(fig_anomaly, f"{name}_anomalies_time_series.png")
        save_figure(fig_field, f"{name}_fields_time_series.png")
        
        plt.close(fig_energy)
        plt.close(fig_anomaly)
        plt.close(fig_field)
        
    if interactive:
        matplotlib.use(original_backend)  # Restore the original backend after processing