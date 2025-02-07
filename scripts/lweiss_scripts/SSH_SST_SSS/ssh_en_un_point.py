import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

data = '/lus/work/CT1/c1601279/lweiss/CROCO/RUN/SWIOSE/OUTPUT/'
work = '/lus/work/CT1/c1601279/lweiss/CROCO/'
path_fig = '/lus/home/CT1/c1601279/lweiss/PYTHON/FIGURES/'

simu = 'run_swio_stogen_1mth'
simu = 'run_swio_stogen_1yr_diff_100p_10jo1ud12'

# Définir les indices du point spécifique
depth_index = 1
lon_index_a = 83
lat_index_a = 171

lon_index_b = 290
lat_index_b = 151

lon_index_c = 231
lat_index_c = 328

# Charger les données
ds1 = xr.open_dataset(data + simu + '_0001/001swiose_his.nc')  # , engine='h5netcdf')
# ssh1 = ds1['zeta'][:, depth_index, lat_index, lon_index]  # Température au point spécifié
time1 = ds1['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
lon_a, lat_a = ds1['lon_rho'][lat_index_a, lon_index_a], ds1['lat_rho'][lat_index_a, lon_index_a]
lon_b, lat_b = ds1['lon_rho'][lat_index_b, lon_index_b], ds1['lat_rho'][lat_index_b, lon_index_b]
lon_c, lat_c = ds1['lon_rho'][lat_index_c, lon_index_c], ds1['lat_rho'][lat_index_c, lon_index_c]
ds1.close()

ds2 = xr.open_dataset(data + simu + '_0002/002swiose_his.nc')  # , engine='h5netcdf')
# ssh2 = ds2['zeta'][:, depth_index, lat_index, lon_index]  # Température au point spécifié
time2 = ds2['time'][:]  # Supposons que 'time' est la variable de temps dans le dataset
ds2.close()

# Convertir les temps en format compréhensible si nécessaire (par exemple, en utilisant pandas)
# time = pd.to_datetime(time, unit='s')  # Ajustez le format de temps si nécessaire

# Créer la figure
fig, ax = plt.subplots(3, 1, figsize=(5, 10), tight_layout=True)

days = np.arange(0, len(time1))
D = np.arange(1, int(len(time1)/24) + 2, 30)
mths = ['J¹⁷', 'F¹⁷', 'M¹⁷', 'A¹⁷', 'M¹⁷', 'J¹⁷', 'J¹⁷', 'A¹⁷'] #, 'S¹⁷', 'O¹⁷', 'N¹⁷', 'D¹⁷'] 

# Tracer la température contre le temps
# ax.plot(time, temp, marker='o', linestyle='-', color='b')
ax[0].plot(days, ds1['zeta'][:, lat_index_a, lon_index_a] - ds2['zeta'][:, lat_index_a, lon_index_a], color='k', linewidth=0.8,linestyle='-') # , label="member 1")
# ax[0].plot(days, ds2['zeta'][:, lat_index_a, lon_index_a], color='r', linewidth=0.8,linestyle='-', label="member 2")
ax[0].set_ylabel('Sea Surface Height Bias [m]')
ax[0].set_xticks(days[0:len(time1):48])
ax[0].set_xticklabels(D, rotation=0)
# ax[0].set_xticklabels(mths, rotation=0)
ax[0].grid(True, linestyle='--', alpha=0.3)
ax[0].set_xlim(days[0], days[-1] + 1)
# ax[0].legend(loc=0, ncol=1, fontsize=10)
ax[0].set_title(f'SSH bias P1, lat {lat_a:.2f}, lon {lon_a:.2f}')
#ax.set_xlabel('Time')

ax[1].plot(days, ds1['zeta'][:, lat_index_b, lon_index_b] - ds2['zeta'][:, lat_index_b, lon_index_b], color='k', linewidth=0.8,linestyle='-') #, label="member 1")
# ax[1].plot(days, ds2['zeta'][:, lat_index_b, lon_index_b], color='r', linewidth=0.8,linestyle='-', label="member 2")
ax[1].set_ylabel('Sea Surface Height Bias [m]')
ax[1].set_xticks(days[0:len(time1):48])
ax[1].set_xticklabels(D, rotation=0)
ax[1].grid(True, linestyle='--', alpha=0.3)
ax[1].set_xlim(days[0], days[-1] + 1)
ax[1].set_title(f'SSH bias P2, lat {lat_b:.2f}, lon {lon_b:.2f}')

ax[2].plot(days, ds1['zeta'][:, lat_index_c, lon_index_c] - ds2['zeta'][:, lat_index_c, lon_index_c], color='k', linewidth=0.8,linestyle='-') #, label="member 1")
# ax[2].plot(days, ds2['zeta'][:, lat_index_c, lon_index_c], color='r', linewidth=0.8,linestyle='-', label="member 2")
ax[2].set_ylabel('Sea Surface Height Bias [m]')
ax[2].set_xticks(days[0:len(time1):48])
ax[2].set_xticklabels(D, rotation=0)
ax[2].set_xticklabels(mths, rotation=0)i
ax[2].grid(True, linestyle='--', alpha=0.3)
ax[2].set_xlim(days[0], days[-1] + 1)
ax[2].set_title(f'SSH bias P3, lat {lat_c:.2f}, lon {lon_c:.2f}')

plt.tight_layout()
plt.savefig(path_fig + '/ssh_bias_time_serie_grid_cells_' + simu + '.png', dpi=300, bbox_inches='tight')
plt.close()

