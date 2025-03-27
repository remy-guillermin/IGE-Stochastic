import cmocean
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

data_path = '/home/guilremy/IGE-Stochastic/Data/swio_his.nc'

visc_cmap = cmocean.cm.amp
diff_cmap = cmocean.cm.tempo

ds = xr.open_dataset(data_path)
AKv = ds.AKv
AKt = ds.AKt
AKs = ds.AKs
s_rho = ds.s_rho
Cs_rho = ds.Cs_rho
hc = ds.hc.values
ds.close()

date_index = -1
date = np.datetime_as_string(AKv.time.data[date_index], 'h')

s_rho_index = -1

dates_index = (-36, -1)
dates = np.datetime_as_string(AKv.time.data[dates_index[0]: dates_index[1]], 'D')

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
fig.suptitle(rf'Vertical Coefficients PDF - {dates[0]} to {dates[-1]} - $\sigma=${s_rho[s_rho_index]:.2f}')

# Temperature Diffusion ($AKt$)
ax = axes[0]
ax.set_title('Temperature Diffusion ($AKt$)')

data = AKt[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
ax.hist(data, bins, color='blue', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKt$')

# Viscosity ($AKv$)
ax = axes[1]
ax.set_title('Viscosity ($AKv$)')

data = AKv[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
ax.hist(data, bins, color='green', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKv$')

# Salinity Diffusion ($AKs$)
ax = axes[2]
ax.set_title('Salinity Diffusion ($AKs$)')

data = AKs[dates_index[0]: dates_index[1], s_rho_index, :, :]
data = data.where((data != 0), np.nan)

# Flatten the data, drop NaNs, and apply log10
data = data.values.flatten()
data = data[~np.isnan(data)]
data = data[data > np.nanmin(data)]
print(f'AKt: min={np.nanmin(data):.2e}, max={np.nanmax(data):.2e}')

bins = np.logspace(np.log10(np.nanmin(data)), np.log10(np.nanmax(data)), 50)
# Plot the PDF using a histogram
n, _, _ = ax.hist(data, bins, color='red', label='Data')
ax.hist(data, bins, color='black', histtype='step')
ax.set_xlabel(r'$AKs$')


for ax in axes:
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.legend()

axes[0].set_ylabel('Count')
plt.tight_layout()

plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/coefficients_PDF_{len(dates)}days.png', dpi=300, transparent=True)
plt.show()
