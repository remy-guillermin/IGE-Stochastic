from setuptools import setup, find_packages

setup(
    name='croco_plot',
    version='0.0.3',
    packages=find_packages(where='modules'),
    package_dir={'': 'modules'},
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
