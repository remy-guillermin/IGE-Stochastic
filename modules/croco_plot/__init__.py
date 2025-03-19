"""
A Python package designed to facilitate the generation of plots for data analysis in the Croco project. This package provides utilities for visualizing data in various formats such as `numpy` and `xarray`.

Main Features:
- Plot climate and oceanographic data.
- Supports visualization libraries such as `matplotlib`, `cartopy`, and `cmocean`.

Modules:
- utils.py: Contains basic functions, including loading grid or his-type files.
- time_series.py: Contains functions to compute and display time series of turbulent kinetic energy (EKE).
- agulhas.py: Contains functions to display zoomed-in values for the Agulhas Current region.
- dataset.py: Contains functions to manipulate netCDF files.
- spectrum.py: Contains functions to display spectra.
- physical_fields.py: Contains functions to analyze and visualize physical fields such as temperature and salinity.
- dynamic_fields.py: Contains functions to analyze and visualize dynamic fields such as velocity and vorticity.

Installation:
    pip install -e modules/

Dependencies:
    numpy, xarray, matplotlib, cmocean, cmcrameri, cartopy, metpy
"""

__version__ = '0.6.0'

from . import utils
from . import dynamic_fields
from . import physical_fields
from . import time_series
from . import agulhas
from . import dataset
from . import spectrum