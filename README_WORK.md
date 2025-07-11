# $WORKDIR
Here are all the data after post-processing (NaN removal, computing of the MLD and wind stress)

## CROCO

This folder contains `croco` source code that has been used to discover the model.

## CROCO_dev_sto_CD

This folder contains a copy of Lisa's `croco` modified code with a perturbation of the wind stress. I used it to develop the vertical turbulent mixing scheme perturbation.

## CROCO_dev_sto_GLS

This folder contains the `croco` with the vertical turbulent mixing scheme perturbation.

## DATASETS

This folder contains data from the simulation after post-processing et contains only temporal series of the mean value in each zone for the Salinity, Temperature and Free Surface.

There is three folder containing valeus for the simulation, de GLORYS reanalysis and different observation sources.

## grid 

This folder contains differents grid files (Horizontal, Vertical)

## MLD

This folder contains `.nc` files post-processed for the MLD computing. Both the interpolated value and the above level value are saved.


## MLD_old

This folder contains the first iteration of `.nc` files post-processed for the MLD computing. These values are wrong because of a mistake in the processing code that did not take into account the local density variation.

## NaN_CORRECTED

This folder contains temperature, salinity and free surface field with the NaN value added in place of the placeholder value.

## NaN_MERGED

Same as before but this time the file are merged for a three year span. There is also ensemble mean and standard deviation files.

## OBS

This folder contains observations from diverse products of the CMEMS post-processed for our need (removal of unwanted variables and unwanted regions).

## RAW

This folder contains raw observation data.

## REGRIDDED

This folder contains interpolated observation on the same grid as the simulation.

# WINDSTR

This folder contains `.nc` files of the post-processed datas for the computing of the wind stress. There is both directional values and the norm of the wind stress.