# ensembles

In this folder are scripts used to make some diagnostics on the ensembles.

## Scripts
- [all_ensemble_spread_tracers_boxes.py](./all_ensemble_spread_tracers_boxes.py): This script plot time series of the tracers (Spatial mean, Ensemble dispersion, RMSE) for each sub-region.
- [all_ensemble_spread_tracers_points.py](./all_ensemble_spread_tracers_points.py): This script plot time series of the tracers (Spatial mean, Ensemble dispersion) for a single point in each sub-region.
- [ensemble_spread_tracers.py](./ensemble_spread_tracers.py): This script plots the mean value for each member of an ensemble for a given year in each sub-region. It also plots the member deviation the same way and the RMSE for each sub-region.
![](../figures/Equator_temp_deviation.png)
![](../figures/Equator_temp.png)
![](../figures/ssh_rms.png)
- [film_std_tracers.py](./film_std_tracers.py): This script use create a [film](https://youtu.be/RAoVWJQzty0) displaying the evolution of the standard deviation for the tracers for any given ensemble. 
- [histogram_dispersion.py](./histogram_dispersion.py): This script compute the histogram of the dispersion for 4 given time steps : 1 day, 10 days, 90 days and 900 days.
![](../figures/windstr_mld_histograms.png)
- [relative_dispersion.py](./relative_dispersion.py): This scripts plots the dispersion of the simulation for each sub-region as a time series, with a normalization corresponding to the dispersion one day after the beginning of the simulation. It also outputs the doubling times, not interpolated.

<table style="width:100%">
    <tr>
        <th rowspan="2">Zone</th>
        <th rowspan="2">Ensemble</th>
        <th colspan="4" style="text-align: center;">Variable</th>
    </tr>
    <tr>
        <td> Sea Surface Salinity </td>
        <td> Sea Surface Temperature </td>
        <td> Mixing Layer Depth </td>
        <td> Wind Stress </td>
    </tr>
    <tr>
        <th rowspan="3">Equator</th>
        <th> INI </th>
        <td style="text-align: right;"> 6 days </td>
        <td style="text-align: right;"> 8 days </td>
        <td style="text-align: right;"> 23 days </td>
        <td style="text-align: right;"> 32 days </td>
    </tr>
    <tr>
        <th> STR </th>
        <td style="text-align: right;"> 2 days </td>
        <td style="text-align: right;"> 3 days </td>
        <td style="text-align: right;"> 3 days </td>
        <td style="text-align: right;"> 9 days </td>
    </tr>
    <tr>
        <th> GLS </th>
        <td style="text-align: right;"> 4 days </td>
        <td style="text-align: right;"> 6 days </td>
        <td style="text-align: right;"> 7 days </td>
        <td style="text-align: right;"> 5 days </td>
    </tr>
    <tr>
        <th rowspan="3">Islands</th>
        <th> INI </th>
        <td style="text-align: right;"> 12 days </td>
        <td style="text-align: right;"> 21 days </td>
        <td style="text-align: right;"> 33 days </td>
        <td style="text-align: right;"> 29 days </td>
    </tr>
    <tr>
        <th> STR </th>
        <td style="text-align: right;"> 2 days </td>
        <td style="text-align: right;"> 14 days </td>
        <td style="text-align: right;"> 10 days </td>
        <td style="text-align: right;"> 7 days </td>
    </tr>
    <tr>
        <th> GLS </th>
        <td style="text-align: right;"> 5 days </td>
        <td style="text-align: right;"> 24 days </td>
        <td style="text-align: right;"> 11 days </td>
        <td style="text-align: right;"> 7 days </td>
    </tr>
    <tr>
        <th rowspan="3">South Mozambique</th>
        <th> INI </th>
        <td style="text-align: right;"> 6 days </td>
        <td style="text-align: right;"> 29 days </td>
        <td style="text-align: right;"> 36 days </td>
        <td style="text-align: right;"> 29 days </td>
    </tr>
    <tr>
        <th> STR </th>
        <td style="text-align: right;"> 4 days </td>
        <td style="text-align: right;"> 2 days </td>
        <td style="text-align: right;"> 4 days </td>
        <td style="text-align: right;"> 15 days </td>
    </tr>
    <tr>
        <th> GLS </th>
        <td style="text-align: right;"> 4 days </td>
        <td style="text-align: right;"> 1 days </td>
        <td style="text-align: right;"> 1 days </td>
        <td style="text-align: right;"> 15 days </td>
    </tr>
    <tr>
        <th rowspan="3">Mascarene</th>
        <th> INI </th>
        <td style="text-align: right;"> 5 days </td>
        <td style="text-align: right;"> 10 days </td>
        <td style="text-align: right;"> 23 days </td>
        <td style="text-align: right;"> 33 days </td>
    </tr>
    <tr>
        <th> STR </th>
        <td style="text-align: right;"> 2 days </td>
        <td style="text-align: right;"> 5 days </td>
        <td style="text-align: right;"> 6 days </td>
        <td style="text-align: right;"> 34 days </td>
    </tr>
    <tr>
        <th> GLS </th>
        <td style="text-align: right;"> 8 days </td>
        <td style="text-align: right;"> 25 days </td>
        <td style="text-align: right;"> 6 days </td>
        <td style="text-align: right;"> 34 days </td>
    </tr>
</table>

![](../figures/all_ens_Mascarene_salt_relative_dispersion_self.png)

## FILM

After creating the frames with [film_std.py](film_std.py), one can create the movie associated with 
```bash
ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p {name}.mp4
```