# Libs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import croco_plot as cplot

dataset_path = '/home/guilremy/IGE-Stochastic/Data/datasets'

datasets = cplot.dataset.fetch_datasets(dir=dataset_path)

variables = 'all'
figsize = (14,8)
name = 'Equator'
roll = 9
xlim = (pd.Timestamp('2017-01-01'), pd.Timestamp('2021-08-31'))

variables = ['Fields']
fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)


files = datasets[name]
# Corrected loop
for file in files:
    source = Path(file).parent.name
    print(source)
    f = xr.open_dataset(file)
    time = pd.to_datetime(f.time.data)
    
    if 'GLORYS' in file:
        color = 'red'
    elif 'CROCO' in file:
        color = 'blue'
    elif 'DATA' in file:
        color = 'black'
    else:
        color = 'white'

    if 'Fields' in variables:
        ax = axes[0]
        sss = f.sss.data

        lw = 1.0

        if 'DATA' in file:
            source = f.sss.source
            lw = 2.0
            
        ax.plot(time, sss, color=color, linestyle='--', alpha=0.3)

        rolling_mean = np.convolve(sss, np.ones(roll)/roll, mode='same')
        ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=lw, label=source)
        ax.set_title('Sea Surface Salinity')
        ax.set_ylabel('SSS [$psu$]')
        ax.set_xlabel('Time')

        ax = axes[1]
        sst = f.sst.data

        lw = 1.0

        if 'DATA' in file:
            source = f.sst.source
            lw = 2.0

        ax.plot(time, sst, color=color, linestyle='--', alpha=0.3)
        rolling_mean = np.convolve(sst, np.ones(roll)/roll, mode='same')
        ax.plot(time[int((roll-1)/2):-int((roll-1)/2)], rolling_mean[int((roll-1)/2):-int((roll-1)/2)], color=color, linestyle='-', linewidth=lw, label=source)
        ax.set_title('Sea Surface Temperature')
        ax.set_ylabel('SST [$°C$]')
        ax.set_xlabel('Time')    
    f.close()

fig.suptitle(f'Fields over time for {name}')

if xlim is not None:
    axes[-1].set_xlim(xlim)

for ax in axes:
    ax.legend(loc="upper right")
    offset = (ax.get_ylim()[1] - ax.get_ylim()[0])/20
    
    start = pd.Timestamp('2018-05-28')
    end = pd.Timestamp('2018-10-21')
    mid = start + (end - start)/2
    ax.axvspan(start, end, color='grey', alpha=0.2)
    ax.text(mid, ax.get_ylim()[1] - offset, '2018 \n Monsoon', fontsize=10, verticalalignment='top', horizontalalignment='center')
    ax.text(mid, ax.get_ylim()[0] + offset, '-3.3%', fontsize=10, verticalalignment='bottom', horizontalalignment='center')
    
    start = pd.Timestamp('2019-06-08')
    end = pd.Timestamp('2019-10-15')
    mid = start + (end - start)/2
    ax.axvspan(start, end, color='grey', alpha=0.2)
    ax.text(mid, ax.get_ylim()[1] - offset, '2019 \n Monsoon', fontsize=10, verticalalignment='top', horizontalalignment='center')
    ax.text(mid, ax.get_ylim()[0] + offset, '+25.3%', fontsize=10, verticalalignment='bottom', horizontalalignment='center')
    
    start = pd.Timestamp('2020-06-01')
    end = pd.Timestamp('2020-10-28')
    mid = start + (end - start)/2
    ax.axvspan(start, end, color='grey', alpha=0.2)
    ax.text(mid, ax.get_ylim()[1] - offset, '2020 \n Monsoon', fontsize=10, verticalalignment='top', horizontalalignment='center')
    ax.text(mid, ax.get_ylim()[0] + offset, '+17.4%', fontsize=10, verticalalignment='bottom', horizontalalignment='center')

fig.tight_layout()

plt.savefig(f'/home/guilremy/IGE-Stochastic/figures/{name}_fields_time_series_monsoon.png', dpi=300, transparent=True)
plt.show()
