import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle

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
lon = g['lon_rho'][:, :] # Longitude
lat = g['lat_rho'][:, :] # Latitude
angle = g['angle'][:, :] # Deformation
msk = g['mask_rho'][:, :] # Mask
msk_inv = np.where(msk == 0, msk, np.nan)
g.close()
print('Grid file closed')

# Process each dataset separately
def process_dataset(d):
    zeta = d['zeta'][:, :, :]
    time = d['time'][:] # Date
    return zeta, time

data_files = [data + simu + 'swio_avg_2017.nc', data + simu + 'swio_avg_2018.nc']
sla_results = {name: [] for name in names}
time_results = []

for data_file in data_files:
    print('Opening data file:', data_file)
    d = xr.open_dataset(data_file)
    print('Data file opened')

    zeta, time = process_dataset(d)
    d.close()
    print('Data file closed')

    # Moyenne annuelle
    zeta_yr = np.nanmean(zeta, axis=0)

    SLA = np.array(zeta - zeta_yr)

    for (lon1, lon2, lat1, lat2), name, color in zip(boxes, names, colors):
        box_mask = np.array((lon >= lon1) & (lon <= lon2) & (lat >= lat1) & (lat <= lat2))
        SLA_box = SLA[:,box_mask]

        # Calculate the spatial mean of EKE_box over time
        SLA_box_mean = np.nanmean(SLA_box, axis=1)
        sla_results[name].append(SLA_box_mean)

    time_results.append(time.data)

# Concatenate results
for name in names:
    sla_results[name] = np.concatenate(sla_results[name])
time_results = np.concatenate(time_results)

fig, axes = plt.subplots(len(names), 1, figsize=(12, 2 * len(names)), sharex=True)

for ax, (name, color) in zip(axes, zip(names, colors)):
    ax.plot(time_results, sla_results[name], color=color, linestyle='--', linewidth=1)
    # Apply centered rolling mean with window of 9 indices (4 before and 4 after)
    sla_rolling_mean = np.convolve(sla_results[name], np.ones(15)/15, mode='same')
    ax.plot(time_results, sla_rolling_mean, color=color, linestyle='-', linewidth=1.5, alpha=0.8)
    ax.set_ylabel('SLA [m]')
    ax.set_title(name)

axes[-1].set_xlabel('Time')
fig.suptitle('Sea Level Anomaly Over Time for Different Boxes')
plt.show()
