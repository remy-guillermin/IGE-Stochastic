from setuptools import setup, find_packages

setup(
    name='croco_plot',
    version='0.0.3',
    packages=find_packages(where='croco_plot'),
    package_dir={'': 'croco_plot'},
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
