# MLD

In this folder are all the scripts related to the mixing layer depth. We found computing scripts as well as plotting, etc.

## Scripts
- [save_mld.py](./save_mld.py): This script computes the Mixing Layer Depth (MLD) using a simple masking method.
- [save_mld_interp.py](./save_mld_interp.py): This script is the improved version for computing the MLD. It takes into account the local density variation and use a linear interpolation to get rid of the vertical levels limitation.
- [comp_mld.py](./comp_mld.py): This script plots slices of the MLD. It also plots the bathymetry for this slice.
![](../figures/slices_mld_comp_2017-09-10.png)
![](../figures/slices_mld_comp_sup_2017-09-10.png)
![](../figures/slices_bathy.png)
- [comp_mld_slice.py](./comp_mld_slice.py): This script plots maps of the ensemble mean and standard deviation of the MLD for the three ensembles.
![](../figures/mld_mean_comp_2017-12-15.png)
![](../figures/mld_std_comp_2017-12-15.png)
- [compute_mld.py](./compute_mld.py): This script plots a map and slices of the MLD. It also plots the mask used to compute the MLD.
![](../figures/slices_mask_mld_2017-12-31.png)
- [density_slice.py](./density_slice.py): This script plots a slice of the density.
![](../figures/density_slice.png)
![](../figures/mld_2020-01-01.png)
![](../figures/slices_mld_2017-12-25.png)
- [mld_dispersion.py](./mld_dispersion.py): This script plots a time series of the dispersion of the MLD.
- [mld_hist_from_files.py](./mld_hist_from_files.py): This script plots the histograms of the MLD of the sub-regions.
![](../figures/hist_ensemble_2017-07-21.png)
- [plot_mld.py](./plot_mld.py): This script plots a map of the interpolated MLD.
![](../figures/mld_interp_2020-01-01.png)
- [plot_mld_from_files.py](./plot_mld_from_files.py): This script plots the MLD from the files in which it is saved.
- [plot_mld_points.py](./plot_mld_points.py): This script plots a superposition of the mean MLD for each ensemble in a given point.
![](../figures/points_Islands_mld_comp_sup_time_series.png)
- [plot_mld_slices_from_files.py](./plot_mld_slices_from_files.py): This script plots a comparison between slices and ensemble of the MLD from the files in which they are saved.
![](../figures/slices_mld_comp_2017-09-10.png)
- [plot_time_series_from_files.py](./plot_time_series_from_files.py): This script plots the time series of the MLD from the files in which they are saved.


