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
from .utils import load_grid, load_data, load_GLORYS, compute_dx_dy

def create_TS_from_CROCO(
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
        lon, lat, _, _, _, _, _, _ = load_grid(grid_path)
    else:
        lon, lat, _, _, _, _, _, _ = load_grid()
        
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
        
    time = pd.to_datetime(time.data)
    
    for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
        results = {var: [] for var in variables}
        
        for var in ['mke', 'eke', 'sla', 'ssa', 'sta', 'ssh', 'sss', 'sst']:
            if var in variables:
                for year in np.unique(time.year):
                    year_mask = time.year == year
                    if var == 'mke':
                        var_data = 1 / 2 * (u[year_mask,:-1,:] ** 2 + v[year_mask,:,:-1] ** 2 + w[year_mask,:-1,:-1] ** 2)
                    elif var == 'eke':
                        ut = (u[year_mask] - np.nanmean(u[year_mask], axis=0))
                        vt = (v[year_mask] - np.nanmean(v[year_mask], axis=0))
                        wt = (w[year_mask] - np.nanmean(w[year_mask], axis=0))
                        var_data = 1 / 2 * (ut[:,:-1,:] ** 2 + vt[:,:,:-1] ** 2 + wt[:,:-1,:-1] ** 2)
                        del ut, vt, wt
                    elif var == 'sla':
                        var_data = (zeta[year_mask] - np.nanmean(zeta[year_mask], axis=0))
                    elif var == 'ssa':
                        var_data = (salt[year_mask] - np.nanmean(salt[year_mask], axis=0))
                    elif var == 'sta':
                        var_data = (temp[year_mask] - np.nanmean(temp[year_mask], axis=0))
                    elif var == 'ssh':
                        var_data = zeta[year_mask]
                    elif var == 'sss':
                        var_data = salt[year_mask]
                    elif var == 'sst':
                        var_data = temp[year_mask]
                    
                    box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                    if var in ['mke', 'eke']:
                        box_mask = box_mask[:-1,:-1]
                    
                    box_data = var_data[:,box_mask]
                    box_sum = np.nanmean(box_data, axis=1)
                    results[var].append(box_sum)
                    del var_data, box_data
        
        results = {var: np.concatenate(results[var]) for var in results}
        
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
                    
def create_TS_from_GLORYS(
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
            
def create_spectrum_from_CROCO(
    data_path:str,
    grid_path:str=None,
    num:int=200,
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
):  
    output_dir = f'/lus/work/CT1/c1601279/rguillermin/datasets/SWIO'
    os.makedirs(output_dir, exist_ok=True)
    # Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, _, _, _, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, _, _, _, _ = load_grid()
        
    # Compute effective grid spacing
    dx_eff = np.mean(1 / pm.data)  # Mean dx
    dy_eff = np.mean(1 / pn.data)  # Mean dy
        
    if isinstance(data_path, str):
        data_path = [data_path]
    
    fill_value = 9.96921e+36
    
        
    num_times_tot = 0
    
    results = {name: [] for name in names}
    
    mean_spectrum = np.zeros(num)  # Initialize mean spectrum
    
    for path in data_path:
        u, v = load_data(path, ('u', 'v'))
        u = u[:,-1,:,:].where((u[:,-1,:,:] != fill_value), np.nan)
        v = v[:,-1,:,:].where((v[:,-1,:,:] != fill_value), np.nan)
        print("NaN values added")
        u = u.data
        v = v.data
        print("Data converted to NumPy arrays")

        K_bins_box = {name: [] for name in names}
        K_box = {name: [] for name in names}
        
        for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            true_indices = np.nonzero(box_mask)
            num_true_points = np.sum(box_mask)

            # Get the bounding box of the True region
            if num_true_points > 0:
                min_row, max_row = np.min(true_indices[0]), np.max(true_indices[0])
                min_col, max_col = np.min(true_indices[1]), np.max(true_indices[1])
                mask_shape = (int(max_row - min_row + 1), int(max_col - min_col + 1))
            # Define wavenumber grid (adjusted for physical spacing)
            ny, nx = mask_shape
            kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx_eff))  # Convert 1/meters
            ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy_eff))
            # Compute isotropic wavenumber k
            KX, KY = np.meshgrid(kx, ky)
            K = np.sqrt(KX**2 + KY**2)
            K_box[name] = K
            k_bins = np.linspace(0, np.max(K), num=int(num/2))
            K_bins_box[name] = k_bins
   
        # Define wavenumber grid (adjusted for physical spacing)
        ny, nx = u[0,:-1,].shape
        kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx_eff))  # Convert 1/meters
        ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy_eff))

        # Compute isotropic wavenumber k
        KX, KY = np.meshgrid(kx, ky)
        K = np.sqrt(KX**2 + KY**2)  # Scalar wavenumber

        # Define wavenumber bins for spectrum calculation
        k_bins = np.linspace(0, np.max(K), num=num)  
        
        # Initialize accumulation for mean spectrum
        num_times = u.shape[0]  # Number of time steps
        spectrum_sum = np.zeros_like(k_bins)  # Initialize sum of spectra

        for t in range(num_times):
            if t % 50 == 0:
                print(f"Processing time step {t+1}/{num_times}")

            # Fill NaNs with zero for FFT
            u_filled = np.nan_to_num(u[t, :, :], nan=0.0)
            v_filled = np.nan_to_num(v[t, :, :], nan=0.0)

            # Compute velocity magnitude squared
            velocity = np.sqrt(u_filled[:-1, :]**2 + v_filled[:, :-1]**2)

            # Apply Hann window to reduce spectral leakage
            window = np.outer(np.hanning(velocity.shape[0]), np.hanning(velocity.shape[1]))
            velocity_windowed = velocity * window
            
            for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                true_indices = np.nonzero(box_mask)
                num_true_points = np.sum(box_mask)
                box_mask =  box_mask[:-1,:-1]

                # Get the bounding box of the True region
                if num_true_points > 0:
                    min_row, max_row = np.min(true_indices[0]), np.max(true_indices[0])
                    min_col, max_col = np.min(true_indices[1]), np.max(true_indices[1])
                    mask_shape = (int(max_row - min_row + 1), int(max_col - min_col + 1))
                
                U_k = np.fft.fft2(velocity_windowed[min_row:max_row,min_col:max_col].reshape(mask_shape))
                U_k_shifted = np.fft.fftshift(U_k)
                power_spectrum = np.abs(U_k_shifted) ** 2 / 1000
                # Bin the spectrum in k-space
                spectrum = np.zeros_like(K_bins_box[name])
                for i in range(len(K_bins_box[name]) - 1):
                    mask_k = (K_box[name] >= K_bins_box[name][i]) & (K_box[name] < K_bins_box[name][i + 1])  # Select wavenumber shell
                    spectrum[i] = np.mean(power_spectrum[mask_k])  # Average over shell
                    
                if type(results[name]) == list:
                    results[name] = spectrum
                else:
                    results[name] += spectrum

            # Compute 2D FFT and shift zero-frequency to the center
            U_k = np.fft.fft2(velocity_windowed)
            U_k_shifted = np.fft.fftshift(U_k)

            # Compute the power spectrum (|FFT|^2)
            power_spectrum = np.abs(U_k_shifted) ** 2 / 1000

            # Bin the spectrum in k-space
            spectrum = np.zeros_like(k_bins)
            for i in range(len(k_bins) - 1):
                mask_k = (K >= k_bins[i]) & (K < k_bins[i + 1])  # Select wavenumber shell
                spectrum[i] = np.mean(power_spectrum[mask_k])  # Average over shell

            # Accumulate for mean spectrum
            spectrum_sum += spectrum
    
        num_times_tot +=  num_times
        # Compute mean spectrum over time
        mean_spectrum += spectrum_sum
        print("Computed mean spectrum over time")
        
        for name in names:
            K_bins_box[name] = K_bins_box[name] * 1000  # Convert to m^-1
        k_bins = k_bins * 1000  # Convert to m^-1
    
    mean_spectrum = mean_spectrum / num_times
    for name in names:
        results[name] = np.array(results[name]) / num_times
    for name in names:
        ds = xr.Dataset({
            "mean_spectrum": (["k"], results[name])},
            coords={
                "k": K_bins_box[name]
            },
            attrs={
                "description": "Mean spectrum over time",
            }
        )
        file_path = os.path.join(output_dir, f"box_data_energy_{name}.nc")
        ds.to_netcdf(file_path)
        print(f'File {file_path} saved')
        
    ds = xr.Dataset({
            "mean_spectrum": (["k"], mean_spectrum)},
            coords={
                "k": k_bins
            },
            attrs={
                "description": "Mean spectrum over time",
            }
        )
    file_path = os.path.join(output_dir, f"global_energy.nc")
    ds.to_netcdf(file_path)
    print(f'File {file_path} saved')
    
def create_spectrum_from_GLORYS(
    data_path:str,
    num:int=200,
    boxes:tuple=[(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)], 
    names:tuple=['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
):  
    output_dir = f'/lus/work/CT1/c1601279/rguillermin/datasets/GLORYS'
    os.makedirs(output_dir, exist_ok=True)
        
    if isinstance(data_path, str):
        data_path = [data_path]
        
    num_times_tot = 0
    
    results = {name: [] for name in names}
    
    mean_spectrum = np.zeros(num)
    
    for path in data_path:
        _, _, _, u, v, lon, lat, _, _ = load_GLORYS(path)
        dx_eff, dy_eff = compute_dx_dy(lat, lon)
        dx_eff = np.mean(dx_eff)
        dy_eff = np.mean(dy_eff)
        print("NaN values added")
        u = u[:,0,:,:].data
        v = v[:,0,:,:].data
        print("Data converted to NumPy arrays")
        
        K_bins_box = {name: [] for name in names}
        K_box = {name: [] for name in names}
        
        for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
            box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
            true_indices = np.nonzero(box_mask)
            num_true_points = np.sum(box_mask)

            # Get the bounding box of the True region
            if num_true_points > 0:
                min_row, max_row = np.min(true_indices[0]), np.max(true_indices[0])
                min_col, max_col = np.min(true_indices[1]), np.max(true_indices[1])
                mask_shape = (int(max_row - min_row + 1), int(max_col - min_col + 1))
            # Define wavenumber grid (adjusted for physical spacing)
            ny, nx = mask_shape
            kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx_eff))  # Convert 1/meters
            ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy_eff))
            # Compute isotropic wavenumber k
            KX, KY = np.meshgrid(kx, ky)
            K = np.sqrt(KX**2 + KY**2)
            K_box[name] = K
            k_bins = np.linspace(0, np.max(K), num=int(num/2))
            K_bins_box[name] = k_bins

        # Define wavenumber grid (adjusted for physical spacing)
        ny, nx = u[0,:,:].shape
        kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx_eff))  # Convert 1/meters
        ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy_eff))

        # Compute isotropic wavenumber k
        KX, KY = np.meshgrid(kx, ky)
        K = np.sqrt(KX**2 + KY**2)  # Scalar wavenumber
        
        # Define wavenumber bins for spectrum calculation
        k_bins = np.linspace(0, np.max(K), num=num)  
        
        # Initialize accumulation for mean spectrum
        num_times = u.shape[0]  # Number of time steps
        spectrum_sum = np.zeros_like(k_bins)  # Initialize sum of spectra

        for t in range(num_times):
            if t % 50 == 0:
                print(f"Processing time step {t+1}/{num_times}")

            # Fill NaNs with zero for FFT
            u_filled = np.nan_to_num(u[t, :, :], nan=0.0)
            v_filled = np.nan_to_num(v[t, :, :], nan=0.0)

            # Compute velocity magnitude squared
            velocity = np.sqrt(u_filled[:, :]**2 + v_filled[:, :]**2)

            # Apply Hann window to reduce spectral leakage
            window = np.outer(np.hanning(velocity.shape[0]), np.hanning(velocity.shape[1]))
            velocity_windowed = velocity * window
            
            for (lon1, lon2, lat1, lat2), name in zip(boxes, names):
                box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
                true_indices = np.nonzero(box_mask)
                num_true_points = np.sum(box_mask)

                # Get the bounding box of the True region
                if num_true_points > 0:
                    min_row, max_row = np.min(true_indices[0]), np.max(true_indices[0])
                    min_col, max_col = np.min(true_indices[1]), np.max(true_indices[1])
                    mask_shape = (int(max_row - min_row + 1), int(max_col - min_col + 1))
                U_k = np.fft.fft2(velocity_windowed[box_mask].reshape(mask_shape))
                U_k_shifted = np.fft.fftshift(U_k)
                power_spectrum = np.abs(U_k_shifted) ** 2 / 1000
                # Bin the spectrum in k-space
                spectrum = np.zeros_like(K_bins_box[name])
                for i in range(len(K_bins_box[name]) - 1):
                    mask_k = (K_box[name] >= K_bins_box[name][i]) & (K_box[name] < K_bins_box[name][i + 1])  # Select wavenumber shell
                    spectrum[i] = np.mean(power_spectrum[mask_k])  # Average over shell
                    
                if type(results[name]) == list:
                    results[name] = spectrum
                else:
                    results[name] += spectrum

            # Compute 2D FFT and shift zero-frequency to the center
            U_k = np.fft.fft2(velocity_windowed)
            U_k_shifted = np.fft.fftshift(U_k)

            # Compute the power spectrum (|FFT|^2)
            power_spectrum = np.abs(U_k_shifted) ** 2 / 1000

            # Bin the spectrum in k-space
            spectrum = np.zeros_like(k_bins)
            for i in range(len(k_bins) - 1):
                mask_k = (K >= k_bins[i]) & (K < k_bins[i + 1])  # Select wavenumber shell
                spectrum[i] = np.mean(power_spectrum[mask_k])  # Average over shell

            # Accumulate for mean spectrum
            spectrum_sum += spectrum

        num_times_tot +=  num_times
        # Compute mean spectrum over time
        mean_spectrum += spectrum_sum
        print("Computed mean spectrum over time")
        
        for name in names:
            K_bins_box[name] = K_bins_box[name] * 1000  # Convert to m^-1
        k_bins = k_bins * 1000  # Convert to m^-1
    
    mean_spectrum = mean_spectrum / num_times
    for name in names:
        results[name] = np.array(results[name]) / num_times
    for name in names:
        ds = xr.Dataset({
            "mean_spectrum": (["k"], results[name])},
            coords={
                "k": K_bins_box[name]
            },
            attrs={
                "description": "Mean spectrum over time",
            }
        )
        file_path = os.path.join(output_dir, f"box_data_energy_{name}.nc")
        ds.to_netcdf(file_path)
        print(f'File {file_path} saved')
        
    ds = xr.Dataset({
            "mean_spectrum": (["k"], mean_spectrum)},
            coords={
                "k": k_bins
            },
            attrs={
                "description": "Mean spectrum over time",
            }
        )
    file_path = os.path.join(output_dir, f"global_energy.nc")
    ds.to_netcdf(file_path)
    print(f'File {file_path} saved')

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
    
def fetch_datasets(
    file_pattern: str  = '*.nc'  
) -> dict:
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
    for file in current_directory.rglob(file_pattern):
        # Extract the name part (e.g., 'Equator' from 'box_data_Equator.nc')
        name = file.stem.replace('box_data_', '')  # Remove 'box_data_' from the filename stem

        if name in file_groups:
            file_groups[name].append(str(file))
        else:
            file_groups[name] = [str(file)]
            
    return file_groups