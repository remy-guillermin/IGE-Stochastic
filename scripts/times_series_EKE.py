import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs

# Load data (same as your existing code)
data = '/lus/store/CT1/c1601279/lweiss/RUN_CROCO/'
simu = 'run_swio2_deter_2017_2023_complet/'
grid = '/lus/store/CT1/c1601279/lweiss/GRID/croco_grid_swio2.nc'

# Define zones
boxes = [(48, 60, -4, 3), (41, 47, -15, -8), (36.5, 42.5, -28, -19), (52, 60, -24, -16)]
names = ['Equator', 'Mayotte-Comores', 'South-Moz', 'Mascarene']
colors = ['saddlebrown', 'darkorchid', 'navy', 'teal']

print('Opening grid file:', grid)
g = xr.open_dataset(grid)
h = g['h'][:, :] # Bathymetry
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()
print('Grid file closed')


data_files = [data + simu + 'swio_avg_2017.nc']
eke_results = {name: [] for name in names}
vertical_eke_results = {name: [] for name in names}
time_results = []

for data_file in data_files:
    print('Opening data file:', data_file)
    d = xr.open_dataset(data_file)
    print('Data file opened')
    u = d['u'][:, :, :, :] # Vitesse surface u
    v = d['v'][:, :, :, :] # Vitesse surface v
    w = d['w'][:, :, :, :] # Vitesse surface w
    time = d['time'][:] # Date
    s_rho = d['s_rho'][:] # S coordinate at RHO-points
    d.close()
    print('Data file closed')
    
    depth = h * s_rho # Profondeur
    depth = np.transpose(depth.data, (2, 0, 1)) # Transpose depth to match u, v, w shape
    
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

    print('Transforming data to geographical coordinates')
    # Transformation des composantes de vent (grille déformée -> grille géographique) pour chaque time index
    angle_expand = angle[:,:].data.reshape(1, angle.shape[0], angle.shape[1], 1)
    
    ut_geo = ut[:,:-1,:,:].data * np.cos(angle_expand[:,:-1,:-1:,]) - vt[:,:,:-1,:].data * np.sin(angle_expand[:,:-1,:-1,:])
    vt_geo = ut[:,:-1,:,:].data * np.sin(angle_expand[:,:-1,:-1:,]) + vt[:,:,:-1,:].data * np.cos(angle_expand[:,:-1,:-1,:])
    wt_geo = wt[:,:-1,:-1,:].data

    print('Calculating EKE')
    EKE = 1 / 2 * (ut_geo ** 2 + vt_geo ** 2 + wt_geo ** 2)
    print('shape ', EKE.shape, 'max EKE: ', np.round(float(np.nanmax(EKE)), 3))
    
    print('Integrating EKE over depth')
    depth_expand = depth[:,:-1,:-1].reshape(EKE.shape[0], EKE.shape[1], EKE.shape[2], 1)
    EKE_weighted = EKE * depth_expand
    EKE_integrated = np.sum(EKE_weighted, axis=0)

    print('Extracting EKE for each box')
    for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
        box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))[:-1,:-1]
        EKE_box = EKE[:,box_mask,:]
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

fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

for ax, (name, color) in zip(axes, zip(names, colors)):
    ax.semilogy(time_results, eke_results[name][:,-1], label=name, color=color)
    ax.set_ylabel('EKE [m²/s²]')
    ax.legend()
    ax.set_ylim(1e2,5e3)

axes[-1].set_xlabel('Time')
fig.suptitle('EKE Over Time for Different Boxes')
plt.show()
