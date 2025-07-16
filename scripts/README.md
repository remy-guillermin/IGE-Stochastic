# Scripts
In this folder are the scripts used during the internship.

## Important folders

- [comparison](comparison/README.md): In this folder are two scripts used to plot comparison between ensemble or with the observation.
- [copernicus](copernicus/README.md): In this folder are scripts used to work with the copernicus marine python API to download data.
- [data](data/README.md): In this folder are scripts used to manipulate the large data files.
- [ensembles](ensembles/README.md): In this folder are scripts used to do some ensemble computations on the data, such as mean, spread or histograms.
- [grid](grid/README.md): In this folder are scripts used to work with the grid. We have bathymetry plot, vertical level computation and regridding.
- [MLD](MLD/README.md): In this folder are scripts used to manipulate the Mixing Layer Depth (MLD), mainly computing and plotting.
- [time_series](time_series/README.md): In this folder are scripts used to plot time series.

## Important scripts

- [fetch_CMEMS.py](copernicus/fetch_CMEMS.py): This script is used to fetch and download data from the CMEMS database.
- [regrid.py](grid/regrid.py): This script is used to regrid any datasets onto the simulation grid. It can either be Copernicus or GLORYS data.
- [compute_wind_str.py](data/compute_wind_str.py): This script is used to compute the wind stress and save it to a different `.nc` file.
- [save_mld_inter.py](MLD/save_mld_interp.py): This script is used to compute the MLD both using a linear interpolation and simply a threshold on the depth levels, and save it to a different `.nc` file.
- [dispersion.py](comparison/dispersion.py): This script is used to plot to ensemble dispersion (standard deviation) of a chosen variable (Surface fields, MLD or Wind stress) for a given date.
- [MLD/plot_time_series_from_files.py](scripts/MLD/plot_time_series_from_files.py): This script get the data saved in the files in the WORK directory and plot the spatially mean value for each boxes we chose as a time series.