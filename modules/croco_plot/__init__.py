"""
Un package Python créé pour faciliter la génération de graphiques dans le cadre de l'analyse des données du projet Croco. Ce package contient des utilitaires pour visualiser des données sous différents formats tels que `numpy` et `xarray`.

Fonctionnalités principales :
- Tracer des graphiques de données climatiques et océanographiques.
- Supporte les bibliothèques `matplotlib`, `cartopy` et `cmocean` pour les visualisations.

Modules :
- utils.py : Contient des fonctions de base notamment pour charger des fichiers types grid ou his.
- plot.py : Contient des fonctions pour l'affichage des données CROCO.
- surf_plot.py : Contient des fonctions pour l'affichage des données de surface CROCO.
- time_series.py : Contient des fonctions pour calculer et afficher les séries temporelles de l'énergie cinétique turbulente (EKE).

Installation :
    pip install -e modules/

Dépendances :
    numpy, xarray, matplotlib, cmocean, cmcrameri, cartopy, metpy
"""

__version__ = '0.0.3'

from . import utils
from . import plot
from . import surf_plot
from . import time_series
