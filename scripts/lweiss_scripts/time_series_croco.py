# -------LIB------------------
import numpy as np
import xarray as xr
import sys
import pandas as pd
# -----------------------------

def print_progress(i, max, post_text):
    n_bar = 10  # size of progress bar
    j = i / max
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'*' * int(n_bar * j ):{n_bar}s}] {int(100 * j)}% {post_text}")
    sys.stdout.flush()

lon_min = 33
lon_max = 61 
lat_min = -31
lat_max = 4

    # {'lon_min': 29.3, 'lon_max': 65.1, 'lat_min': -32.9, 'lat_max': 5.1},  # SWIO

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_0003_2017'

###############################################################################
#                                 GRID                                        #
###############################################################################
g = xr.open_dataset(data + simu + '/swiose_grid.nc')
lon_rho = g['lon_rho'][:, :]
lat_rho = g['lat_rho'][:, :]
h = g['h'][:, :]
msk_rho = g['mask_rho'][:, :]
msk_inv = np.where(msk_rho == 0, msk_rho, np.nan)
g.close()
###############################################################################
#                                 SST SSS SLA CROCO                           #
###############################################################################
ds = xr.open_dataset(data + simu + '/swiose_avg.nc')

# Create a mask for the latitude and longitude range
mask = ((ds.lat_rho >= lat_min) & (ds.lat_rho <= lat_max) & (ds.lon_rho >= lon_min) & (ds.lon_rho <= lon_max) & (ds.mask_rho == 1))

# Remplacer les valeurs hors domaine par NaN
fill_value = 9.96921e+36

zeta = ds['zeta'][:, :, :]
zeta = zeta.where((zeta != fill_value) & mask, np.nan)
zeta_mean_yr = zeta.mean(dim='time', skipna=True)
sla_i = zeta - zeta_mean_yr

temp = ds['temp'][:, -1, :, :]
temp = temp.where((temp != fill_value) & mask, np.nan)
salt = ds['salt'][:, -1, :, :]
salt = salt.where((salt != fill_value) & mask, np.nan)
###############################################################################
area = 1 / (g['pm'] * g['pn'])
area = area.where(mask, np.nan)
area2 = xr.DataArray(np.tile(area.expand_dims('time').values, (ds.dims['time'], 1, 1)),
                    dims=['time', 'eta_rho', 'xi_rho'],
                    coords={'time': ds.time,
                            # 's_rho': ds.s_rho,
                            'eta_rho': area.eta_rho,
                            'xi_rho': area.xi_rho})
###############################################################################
sla = sla_i * area2[:, :, :].values
sst = temp * area2[:, :, :].values
sss = salt * area2[:, :, :].values

mean_sla = np.nansum(sla[:,:,:], axis=(1, 2)) / np.nansum(area2[:,:, :].values, axis=(1, 2))
mean_sst = np.nansum(sst[:,:,:], axis=(1, 2)) / np.nansum(area2[:,:, :].values, axis=(1, 2))
mean_sss = np.nansum(sss[:,:,:], axis=(1, 2)) / np.nansum(area2[:,:, :].values, axis=(1, 2))

print(mean_sla.shape, mean_sst.shape, mean_sss.shape)

mean_sla_series = pd.Series(mean_sla)
mean_sla_roll = mean_sla_series.rolling(30).mean()

df = pd.DataFrame({'mean_sst': mean_sst, 'mean_sss': mean_sss, 'mean_sla': mean_sla, 'mean_sla_roll':mean_sla_roll})
print(df, df.shape)
# df = df.dropna()
#print(df, df.shape)
#df = df.drop([59, 60, 426])
#df = df.drop(df.index[426: 446])
#df = df.drop(df.index[0:61])
#print(df, df.shape)

df.to_csv(path_fig + simu + '/time_series/sst_sss_sla_time_series_' + simu + '.csv', index=False)

################################################################################