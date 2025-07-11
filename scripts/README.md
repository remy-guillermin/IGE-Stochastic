# Scripts
In this folder are the scripts used during the internship.

## Important scripts

- [regrid.py](grid/regrid.py): This script is used to regrid any datasets onto the simulation grid. It can either be Copernicus or GLORYS data.
- [compute_wind_str.py](data/compute_wind_str.py): This script is used to compute the wind stress and save it to a different `.nc` file.
- [save_mld_inter.py](MLD/save_mld_interp.py): This script is used to compute the Mixing Layer Depth (MLD) both using a linear interpolation and simply a threshold on the depth levels, and save it to a different `.nc` file.
- 