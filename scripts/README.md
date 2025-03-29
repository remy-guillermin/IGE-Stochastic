# Scripts
Dans ce dossier, je vais regrouper tous les scripts que je vais utiliser à un moment ou à un autre.

# Install CARTOPY offline

To install cartopy offline on ADASTRA (to be able to use cartopy's feature such as land and coastlines), one should manually download all the [Natural Earth shapefiles](https://www.naturalearthdata.com/downloads/) first. In my case I have downloaded LAND, COASTLINES, OCEAN and BOUNDARY LINES in 10m, 50m and 110m as well as MINOR ISLANDS in 10m. The first three and MINOR ISLANDS will need to be extracted in `.cartopy/shapefiles/natural_earth/physical` and the last one in `.cartopy/shapefiles/natural_earth/cultural`.

In the home directory,
```bash
mkdir -p .cartopy/shapefiles/natural_earth/cultural .cartopy/shapefiles/natural_earth/physical
```

One can delete the `html` and `txt` files from the archive when extracted.

Lastly, one need to add the environment variable for cartopy
```bash
echo 'export CARTOPY_DATA_DIR="$HOME/.cartopy"' >> ~/.bashrc
```

One can test if the installation works by running 
```python
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([25, 69, -36, 7])

    land_color = cfeature.COLORS['land']

    minor_islands = cfeature.NaturalEarthFeature(
        category='physical',
        name='minor_islands',
        scale='10m',
        facecolor=land_color,
        edgecolor='black'
    )

    ax.coastlines(resolution='50m')

    ax.add_feature(ccrs.cartopy.feature.LAND, edgecolor='black')
    ax.add_feature(ccrs.cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=0.5)

    ax.add_feature(minor_islands)

    text = AnchoredText('South West Indian Ocean', loc=4, prop={'size': 12}, frameon=True)
    
    ax.add_artist(text)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
```

If an error message mentionning `Downloading: ...50m/cultural/ne_50m_admin_1_states_provinces_lines_shp.zip` or anything like that appears, there was a problem, probably when setting the environment variable.

The following script,
```python
import cartopy
cartopy.config
```

should display
```bash
{
    'pre_existing_data_dir': PosixPath('/lus/home/CT1/c1601279/rguillermin/.cartopy'),
    'data_dir': PosixPath('/lus/home/CT1/c1601279/rguillermin/.local/share/cartopy'),
    'cache_dir': PosixPath('/tmp/cartopy_cache_dir'),
    'repo_data_dir': PosixPath('/lus/home/CT1/c1601279/rguillermin/python_environment/lib/python3.11/site-packages/cartopy/data'),
    'downloaders': {
        (
            'shapefiles',
            'natural_earth'
        ): <cartopy.io.shapereader.NEShpDownloader at 0x7fc091a39b10>,
        (
            'shapefiles',
            'gshhs'
        ): <cartopy.io.shapereader.GSHHSShpDownloader at 0x7fc091a39cd0>
    }
}
```
If `'pre_existing_data_dir': PosixPath('.')`, this means that `CARTOPY_DATA_DIR` was not set successfully.

The structure of the `.cartopy` folder should be
```bash
.cartopy
└── shapefiles
    └── natural_earth
        ├── cultural
        │   ├── ne_10m_admin_0_boundary_lines_land.***
        │   ├── ne_110m_admin_0_boundary_lines_land.***
        │   └── ne_50m_admin_0_boundary_lines_land.***
        └── physical
            ├── ne_10m_coastline.***
            ├── ne_10m_land.***
            ├── ne_10m_minor_islands_coastline.***
            ├── ne_10m_minor_islands.***
            ├── ne_10m_ocean.***
            ├── ne_110m_coastline.***
            ├── ne_110m_land.***
            ├── ne_110m_ocean.***
            ├── ne_50m_coastline.***
            ├── ne_50m_land.***
            └── ne_50m_ocean.***
```

## Notes
### [wind_stress.py](wind_stress.py)
- Lignes 47-48: On doit convertir la grille déformée (qui suit les latitudes/longitudes) en une grille géographique pour l'affichage. Pour cela, on va utiliser $u_{geo} = u \cos{\theta} - v \sin{\theta}$ et $v_{geo} = u \sin{\theta} + v \cos{\theta}$ avec $\theta$ l'angle ... . Comme dit dans l'article [Vogt-Vincent et al.](../bibliography/gmd-16-1163-2023.pdf) section 2.2, on a une grille plus large à l'équateur qu'à la limite sud.
- Ligne 56: `projection=ccrs.PlateCarree()` est une projection équivalente à une carte plate où les latitudes et longitudes sont représentées par une grille régulière.
- Lignes 66 et 69: `transform=ccrs.PlateCarree()` permet d'assurer la bonne interprétation des coordonnées.
