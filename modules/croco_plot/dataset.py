"""
Module plot pour croco_plot.

Ce module contient des fonctions pour la création et l'affichage de fichiers netCDF issues des fichiers CROCO brut.
"""
import os
import xarray as xr
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from .utils import load_grid, load_data, load_GLORYS

def create_from_CROCO(
    data_path:str,
    variables:str = 'all',
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene'], 
    grid_path:str=None
):      
    """
    Create netCDF files from CROCO data for specified variables and geographical boxes.

    Parameters
    ----------
    data_path : str
        Path to the CROCO data file.
    variables : str, optional
        Variables to process, by default 'all'.
    boxes : tuple, optional
        Geographical boxes defined by (lon1, lon2, lat1, lat2), by default specific regions.
    names : tuple, optional
        Names corresponding to the geographical boxes, by default specific names.
    grid_path : str, optional
        Path to the grid data file, by default None.
    """
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
        u, v, w = u[:,-1,:,:].data, v[:,-1,:,:].data, w[:,-1,:,:].data
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
                    var_data = 1 / 2 * (u[:,:-1,:] ** 2 + v[:,:,:-1] ** 2 + w[:,:-1,:-1] ** 2)
                elif var == 'eke':
                    ut = (u - np.nanmean(u, axis=0))
                    vt = (v - np.nanmean(v, axis=0))
                    wt = (w - np.nanmean(w, axis=0))
                    var_data = 1 / 2 * (ut[:,:-1,:] ** 2 + vt[:,:,:-1] ** 2 + wt[:,:-1,:-1] ** 2)
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
            
        
def create_from_GLORYS(
    data_path:str,
    variables:str = 'all',
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
):      
    """
    Create netCDF files from GLORYS data for specified variables and geographical boxes.

    Parameters
    ----------
    data_path : str
        Path to the GLORYS data file.
    variables : str, optional
        Variables to process, by default 'all'.
    boxes : tuple, optional
        Geographical boxes defined by (lon1, lon2, lat1, lat2), by default specific regions.
    names : tuple, optional
        Names corresponding to the geographical boxes, by default specific names.
    """
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
            
def merge(
    file_pattern: str, 
    output_dir: str,
    output_name: str = None
):
    """
    Merge multiple netCDF files matching a pattern into a single file.

    Parameters
    ----------
    file_pattern : str
        Pattern to match the files to concatenate.
    output_dir : str
        Directory to save the concatenated file.
    output_name : str
        Name of the concatenated file.
    """
    os.makedirs(f'/lus/work/CT1/c1601279/rguillermin/{output_dir}', exist_ok=True)
    if output_name is None:
        output_name = file_pattern
    
    file_list = sorted(glob.glob(f'{file_pattern}*.nc'))
    if not file_list:
        print(f"No files found matching {file_pattern}*.nc")
        return

    ds = xr.open_mfdataset(file_list, combine='by_coords')
    output_file = f"/lus/work/CT1/c1601279/rguillermin/{output_dir}/{output_name}.nc"
    
    print(f'{len(file_list)} files matching {output_name} concatenated.')
    ds.to_netcdf(output_file)
    ds.close()
    
    print(f"Concatenation complete: {output_name}")
    
def fetch_datasets() -> dict:
    """
    Fetch all netCDF files in the current directory and group them by their dataset names.

    Returns
    -------
    dict
        A dictionary where keys are dataset names and values are lists of file paths.
    """
    current_directory = Path.cwd()
    file_groups = {}

    # Iterate over all .nc files
    for file in current_directory.rglob("*.nc"):
        # Extract the name part (e.g., 'Equator' from 'box_data_Equator.nc')
        name = file.stem.replace('box_data_', '')  # Remove 'box_data_' from the filename stem

        if name in file_groups:
            file_groups[name].append(str(file))
        else:
            file_groups[name] = [str(file)]
            
    return file_groups