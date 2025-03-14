import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cmcrameri
import cartopy.crs as ccrs
from .utils import load_grid, load_data, save_figure, plot_map


def eddy(
	data_path: str,
	figsize: tuple = (10, 8), 
	grid_path: str = None
):
	# Load grid data
    if grid_path is not None:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid(grid_path)
    else:
        lon, lat, pm, pn, msk, msk_inv, angle, _ = load_grid()
        
    # Load simulation data
    u, v, w = load_data(data_path, ('u', 'v', 'w'))
        
    u = u[:,-1,:,:]
    v = v[:,-1,:,:]
    w = w[:,-1,:,:]  
    
    fill_value = 9.96921e+36
    u = u.where((u != fill_value), np.nan).data
    v = v.where((v != fill_value), np.nan).data
    w = w.where((w != fill_value), np.nan).data
    print("NaN values added")    
    
    angle = angle.data

    u_geo = u[:,:,:-1,:] * np.cos(angle[:-1,:-1]) - v[:,:,:,:-1] * np.sin(angle[:-1,:-1])
    v_geo = u[:,:,:-1,:] * np.sin(angle[:-1,:-1]) + v[:,:,:,:-1] * np.cos(angle[:-1,:-1])
    w_geo = w[:,:,:-1,:-1]    
    
    # Calculate derivatives
    dv_dlon = np.gradient(v_geo, axis=3) * pm[:-1,:-1]
    du_dlat = np.gradient(u_geo, axis=2) * pn[:-1,:-1]
        
    # Calculate vorticity
    vorticity = dv_dlon - du_dlat
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
