# -------LIB------------------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob
import sys
import pandas as pd
import xarray as xr

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

# simu = 'run_swio_0003_2017'
simu = 'run_swio_stogen_1mth_0001'

lon_min = 33
lon_max = 61
lat_min = -31
lat_max = 4

### open netcdf file and variables
### open grid
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon = g['lon_rho'][:, :]
lat = g['lat_rho'][:, :]
msk = g['mask_rho'][:, :]
h = g['h'][:, :]
#sys.exit()
g.close()

msk_inv = np.where(msk == 0, msk, np.nan)

# ds = xr.open_dataset(data + simu + '/swiose_avg.nc')  # , engine='h5netcdf')
ds = xr.open_dataset(data + simu + '/001swiose_his.nc')
u = ds['u'][:, :, :-1, :]  # Sélectionne le premier niveau s_rho
v = ds['v'][:, :, :, :-1]
s_rho = ds['s_rho'][:] # s_rho(s_rho) S-coordinate at RHO-points
Cs_rho = ds['Cs_rho'][:] # Cs_rho(s_rho) S-coordinate stretching curves at RHO-points
hc = ds['hc'].values # S-coordinate parameter, critical depth
# lon = ds['lon_rho'][:, :]
# lat = ds['lat_rho'][:, :]
# msk = ds['mask_rho'][:, :]

mask = ((ds.lat_rho >= lat_min) & (ds.lat_rho <= lat_max) & (ds.lon_rho >= lon_min) & (ds.lon_rho <= lon_max) & (ds.mask_rho == 1))

# Calcul z0 a nonlinear vertical transformation 
def calc_depth(s, Cs, hc, h):
    N = len(s_rho)
    M, L = h.shape
    z0 = np.zeros((N, M, L))
    depth = np.zeros((N, M, L))
    for k in range(N):
        z0[k, :, :] = (hc * s[k] + h * Cs[k]) / (hc + h)
#        depth[k, :, :] = z0[k, :, :] * h ## (hc * s[k] + h * Cs[k])
    return z0

area = 1 / (g['pm'] * g['pn'])
area = area.where(mask, np.nan)
area2 = xr.DataArray(np.tile(area.expand_dims('time').values, (ds.dims['time'], 1, 1)), 
                    dims=['time', 'eta_rho', 'xi_rho'], 
                    coords={'time': ds.time, 
                            # 's_rho': ds.s_rho,
                            'eta_rho': area.eta_rho, 
                            'xi_rho': area.xi_rho})

# area = np.tile(area.values[np.newaxis, :, :], (ds.dims['time'], 1, 1))
# area = area.expand_dims('time', axis=0).broadcast_like(ds['time'])
vol = np.tile(area.values[np.newaxis, :, :], (50, 1, 1)) * calc_depth(s_rho, Cs_rho, hc, h)
vol2 = np.tile(vol[np.newaxis, :, :, :], (u.shape[0], 1, 1, 1))

fill_value = 9.96921e+36
u = u.where(u != fill_value, np.nan)
v = v.where(v != fill_value, np.nan)

ke = 0.5 * (u[:,-1,:,:].values ** 2 + v[:,-1,:,:].values ** 2) * area2[:,:-1, :-1].values
ke_3d = 0.5 * (u[:,:,:,:].values ** 2 + v[:,:,:,:].values ** 2) * vol2[:,:,:-1, :-1]

#ke[np.isinf(ke)] = np.nan
#ke_3d[np.isinf(ke_3d)] = np.nan

mean_ke = np.nansum(ke[:,:,:], axis=(1, 2)) / np.nansum(area2[:,:-1, :-1].values, axis=(1, 2))
mean3d_ke = np.nansum(ke_3d[:,:,:,:], axis=(1, 2, 3)) / np.nansum(vol2[:,:,:-1, :-1], axis=(1, 2, 3))

df = pd.DataFrame({'ke': mean_ke, 'ke_3d': mean3d_ke})

df.to_csv(path_fig + '/ke_time_serie_' + simu + '.csv', index=False)

sys.exit()

###################################################################################################################
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5), tight_layout=True)

df = pd.read_csv(path_fig + 'run_swio_0003_2017/KE/ke_time_serie_run_swio_0003_2017.csv')
df1 = pd.read_csv(path_fig + 'run_swio_stogen_1mth/ke_time_serie_run_swio_stogen_1mth_0001.csv') 
df2 = pd.read_csv(path_fig + 'run_swio_stogen_1mth/ke_time_serie_run_swio_stogen_1mth_0002.csv')

days = np.arange(0, len(df1))
mths = ['J¹⁷', 'F¹⁷', 'M¹⁷', 'A¹⁷', 'M¹⁷', 'J¹⁷', 'J¹⁷', 'A¹⁷', 'S¹⁷', 'O¹⁷', 'N¹⁷', 'D¹⁷']
D = np.arange(1, int(len(df1)/24) + 2, 2) 

# ax[0].plot(days[0:len(df)], df.ke*1000, color='k', linewidth=0.8,linestyle='-', label="croco.SWIO")
ax[0].plot(days[0:len(df1)], df1.ke*1000, color='b', linewidth=0.8,linestyle='-', label="member 1")
ax[0].plot(days[0:len(df2)], df2.ke*1000, color='r', linewidth=0.8,linestyle='-', label="member 2")
ax[0].set_ylabel('Kinetic Energy (10⁻³ m² s⁻²)')
# ax[0].set_xticks(days[15:452:30])
ax[0].set_xticks(days[0:len(df1):48])
# ax[0].set_xticklabels(mths, rotation=0)
ax[0].set_xticklabels(D, rotation=0)
ax[0].grid(True, linestyle='--', alpha=0.3)
ax[0].set_xlim(days[0], days[-1] + 1)
ax[0].legend(loc=0, ncol=1, fontsize=10)
fig.text(0.09, 0.97, "(a) surface", fontsize=12)

# ax[1].plot(days[0:len(df)], df.ke_3d*1000, color='k', linewidth=0.8, linestyle='-', label="croco.SWIO")
ax[1].plot(days[0:len(df1)], df1.ke_3d*1000, color='b', linewidth=0.8, linestyle='-', label="member 1")
ax[1].plot(days[0:len(df2)], df2.ke_3d*1000, color='r', linewidth=0.8, linestyle='-', label="member 2")
# ax[1].set_xticks(days[15:452:30])
ax[1].set_xticks(days[0:len(df1):48])
# ax[1].set_xticklabels(mths, rotation=0)
ax[1].set_xticklabels(D, rotation=0)
ax[1].grid(True, linestyle='--', alpha=0.3)
# ax.set_yticks(np.arange(-2, 3, 1))
# ax.set_ylim(-2.4, 2.4)
ax[1].set_xlim(days[0], days[-1] + 1)
ax[1].legend(loc=0, ncol=1, fontsize=10)
fig.text(0.56, 0.97, "(b) 3 dimensions ", fontsize=12)

# Adjust spacing between subplots
plt.tight_layout()

# Display the plot
plt.savefig(path_fig + '/ke_time_serie_' + simu + '.png', format='png', dpi=300, bbox_inches='tight') #, pad_inches=0.2)
# plt.show()
plt.close()
