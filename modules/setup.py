from setuptools import setup, find_packages

setup(
    name='croco_plot',                # Nom de ton package
    version='0.1',                   # Version du package
    packages=find_packages(),        # Trouver tous les sous-dossiers contenant __init__.py
    install_requires=[               # Dépendances du package
        'numpy',
        'xarray',
        'matplotlib',
        'cmocean',
        'cmcrameri',
        'cartopy',
        'metpy',
    ],
)