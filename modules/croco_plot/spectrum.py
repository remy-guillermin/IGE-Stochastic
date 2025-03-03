"""
Module plot pour croco_plot.

Ce module contient des fonctions pour l'affichage de spectres de données CROCO.
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from .utils import save_figure, get_name_from_filename, get_color_from_filename

def Ekdk(
    datasets: str, 
    figsize: tuple = (10, 5)
):
    if isinstance(datasets, str):
        datasets = [datasets]
        
    # Plot the isotropic energy spectrum
    fig, ax = plt.subplots(1, figsize=figsize, sharex=True)
         
    for file in datasets:
        f = xr.open_dataset(file)
        k_bins = f.k.data
        Ekdk = f.mean_spectrum.data
        f.close()
        
        color = get_color_from_filename(filename=file)
        name = get_name_from_filename(filename=file)
        ax.loglog(k_bins[1:-1], Ekdk[1:-1], color=color, label=name)
        
    def one_over(x):
        """Vectorized 1/x, treating x==0 manually"""
        x = np.array(x, float)
        near_zero = np.isclose(x, 0)
        x[near_zero] = np.inf
        x[~near_zero] = 1 / x[~near_zero]
        return x

    inverse = one_over

    # Define reference slopes for k^-5/3 and k^-3
    k_ref = k_bins[1:-1]  # Avoid zero since log-log plots can't handle k=0
    E_k_5_3 = (k_ref / k_ref[0]) ** (-5/3)  # Normalize at first wavenumber
    E_k_3 = (k_ref / k_ref[0]) ** (-3)  # Normalize at first wavenumber
    
    secax = ax.secondary_xaxis('top', functions=(one_over, inverse))
    secax.set_xlabel(r'$\ell$ [km]')
    
    # Add k^-5/3 and k^-3 reference lines
    ax.loglog(k_ref, E_k_5_3, '--', label=r"$k^{-5/3}$", color="gray", linewidth=1.5)
    ax.loglog(k_ref, E_k_3, '--', label=r"$k^{-3}$", color="red", linewidth=1.5)
    
    ax.set_xlabel(r"$k_{\ell}$ ($km^{-1}$)")
    ax.set_ylabel(r'$\overline{\overline{E(k_{\ell})}}$')
    ax.legend()
    
    fig.suptitle("Mean Isotropic Energy Spectrum")
    plt.tight_layout()
    save_figure(fig, f"Mean_Ekdk.png")
    plt.close(fig)
