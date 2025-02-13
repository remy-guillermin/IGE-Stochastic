from setuptools import setup, find_packages

setup(
    name='croco_plot',
    version='0.0.3',
    author='Remy Guillermin',
    description='A package for plotting CROCO simulation data in 2D and 3D.',
    packages=find_packages(where='.'),
    package_dir={'': '.'},
    install_requires=[
        'numpy',
        'xarray',
        'matplotlib',
        'cmocean',
        'cmcrameri',
        'cartopy',
        'metpy'
    ],
)
