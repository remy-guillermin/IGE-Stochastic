"""
Module plot pour croco_plot.

Ce module contient des fonctions pour la création et l'affichage de fichiers netCDF issues des fichiers CROCO brut.
"""
import os
from pathlib import Path
import xarray as xr
import glob
import numpy as np
import pandas as pd
from .utils import load_grid, load_data, load_GLORYS

def create_from_CROCO(
    data_path:str,
    variables:str = 'all',
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
    grid_path:str=None
):      
    output_dir = f'/lus/work/CT1/c1601279/rguillermin/datasets/{os.path.splitext(data_path)[0]}'
    os.makedirs(output_dir, exist_ok=True)
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['eke', 'mke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']
        
    print(f'Computing {", ".join(var for var in variables)}.')
    
    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, _, _, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, _, _, _ = load_grid()
        
    cell_surface = 1 / (pm * pn).data
    cell_surface[(1-msk).astype(int)] = np.nan
    
    zeta, temp, salt, u, v, w, time, _ = load_data(data_path, ('zeta', 'temp', 'salt', 'u', 'v', 'w', 'time', 's_rho'))

        
    fill_value = 9.96921e+36
    
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
        temp = temp[:,-1,:,:].data
        temp[temp == fill_value] = np.nan
    if 'ssa' in variables or 'sss' in variables:
        print('Converting salinity')
        salt = salt[:,-1,:,:].data
        salt[salt == fill_value] = np.nan
        
    for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
        results = {var: [] for var in variables}
        
        for var in ['mke', 'eke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']:
            if var in variables:
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
                
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                if var in ['mke', 'eke']:
                    box_mask = box_mask[:-1,:-1]
                if var in ['sla', 'ssh', 'sta', 'ssa', 'sss', 'sst']:
                    box_data = var_data[:,box_mask]
                else:
                    box_data = var_data[:,:,box_mask]
                if var in ['mke', 'eke']:
                    box_sum = np.nansum(box_data[:,-1,:], axis=1)
                else:
                    box_sum = np.nanmean(box_data, axis=1)
                results[var].append(box_sum)
                del var_data, box_data
        
        results = {var: np.array(results[var]).squeeze() for var in results}
        
        ds = xr.Dataset({
            var: (["time"], results[var]) for var in results},
            coords={
                "time": time
            },
            attrs={
                "description": "Dataset for specific values in a box",
                "box": name
            },
        )
        file_path = os.path.join(output_dir, f"box_data_{name}.nc")
        ds.to_netcdf(file_path)
        print(f'File {file_path} created')
        del ds, results
            
        
def create_from_GLORYS(
    data_path:str,
    variables:str = 'all',
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
):      
    output_dir = f'/lus/work/CT1/c1601279/rguillermin/datasets/{os.path.splitext(data_path)[0]}'
    os.makedirs(output_dir, exist_ok=True)
    if isinstance(variables, str):
        variables = [variables]
        
    if 'all' in variables:
        variables = ['mke', 'eke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']
        
    print(f'Computing {", ".join(var for var in variables)}.')
        
    salt, temp, zeta, u, v, lon, lat, _, _ = load_GLORYS(data_path)
    time = pd.to_datetime(salt['time'].data)
    u = u[:,0,:,:]
    v = v[:,0,:,:]
    salt = salt[:,0,:,:]
    zeta = zeta
    temp = temp[:,0,:,:]
    
    if 'eke' in variables or 'mke' in variables:
        print('Converting velocities')
        u, v = u.data, v.data
    if 'sla' in variables or 'ssh' in variables:
        print('Converting free surface')
        zeta = zeta.data
    if 'sta' in variables or 'sst' in variables:
        print('Converting temperature')
        temp = temp.data
    if 'ssa' in variables or 'sss' in variables:
        print('Converting salinity')
        salt = salt.data
        
    for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
        results = {var: [] for var in variables if var != 'ke'}
        
        for var in ['mke', 'eke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']:
            if var in variables:
                if var == 'mke':
                    var_data = 1 / 2 * (u ** 2 + v ** 2)
                elif var == 'eke':
                    ut = (u - np.nanmean(u, axis=0))
                    vt = (v - np.nanmean(v, axis=0))
                    var_data = 1 / 2 * (ut ** 2 + vt ** 2 )
                    del ut, vt
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
                
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                if var in ['mke', 'eke']:
                    box_mask = box_mask[:,:]
                if var in ['sla', 'ssh', 'sta', 'ssa', 'sss', 'sst']:
                    box_data = var_data[:,box_mask]
                else:
                    box_data = var_data[:,box_mask]
                    
                box_sum = np.nanmean(box_data, axis=1)
                results[var].append(box_sum)
                del var_data, box_data
        
        results = {var: np.array(results[var]).squeeze() for var in results}
        
        ds = xr.Dataset({
            var: (["time"], results[var]) for var in results},
            coords={
                "time": time
            },
            attrs={
                "description": "Dataset for specific values in a box",
                "box": name
            },
        )
        file_path = os.path.join(output_dir, f"box_data_{name}.nc")
        ds.to_netcdf(file_path)
        print(f'File {file_path} created')
        del ds, results
            
def concatenate(
    file_pattern: str, 
    output_dir: str
):
    os.makedirs(f'/lus/work/CT1/c1601279/rguillermin/{output_dir}', exist_ok=True)
    
    file_list = sorted(glob.glob(f'{file_pattern}*.nc'))
    if not file_list:
        print(f"No files found matching {file_pattern}*.nc")
        return

    ds = xr.open_mfdataset(file_list, combine='by_coords')
    output_file = f"/lus/work/CT1/c1601279/rguillermin/{output_dir}/{file_pattern}.nc"
    
    print(f'{len(file_list)} files matching {file_pattern} concatenated.')
    ds.to_netcdf(output_file)
    ds.close()
    
    print(f"Concatenation complete: {output_file}")