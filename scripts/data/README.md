# data

In this folder are regrouped scripts used to perform various operations on the simulation.

## Scripts

- [compute_wind_str.py](./compute_wind_str.py): This script is used to compute and save the wind stress coefficient for each simulation. It outputs a `.nc` file containing all the three year of data for the whole domain.
- [compute_wind_str_mean.py](./compute_wind_str_mean.py): This script computes the mean and standard deviation of the wind stress and saves it to two new files.
- [remove_NaN.py](./remove_NaN.py): This scripts opens up every average `.nc` file from the simulation, replace the fill_value by NaN and saves those datasets as new files.
- [merge_NaN.py](./merge_NaN.py): This scripts merges all the files created by [remove_NaN.py](./remove_NaN.py) for each member of each simulation.
- [compute_trac_mean.py](./compute_wind_str_mean.py): This script computes the mean and standard deviation of the tracers and saves it to two new files.
- [create_datasets.py](./create_datasets.py): This script create datasets for each chosen sub-region containing the spatially mean value of each tracers as a time series.
- [merge_obs.py](./merge_obs.py): This script merges the observation datasets.
- [merge_swio.py](./merge_swio.py): This script merge the simulation data.