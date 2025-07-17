# Documentation

## Adastra

### Login and configuration

#### Configuration of `./ssh/config`
```bash
Host *
	ServerAliveInterval 30
	ForwardX11 yes

Host ige-ssh
	HostName ige-ssh.u-ga.fr
	User guilremy

Host adastra
	Hostname adastra.cines.fr
	User rguillermin
	ProxyCommand ssh ige-ssh nc %h %p  
```

```bash
ssh adastra
```
We need IGE's password as well as ADASTRA's one.

#### `SSH` configuration to connect to IGE computer with a key.
1. Generate a SSH key:
    ```bash
    ssh-keygen -t rsa -b 4096 -C "remy.guillermin@email.com" -f ~/.ssh/ige_ssh
    ```
2. Copie the SSH key on the IGE computer :
    ```bash
    ssh-copy-id -i ~/.ssh/ige_ssh.pub guilremy@ige-ssh.u-ga.fr
    ```
3. Modify `./ssh/config` :
    ```bash
    Host ige-ssh
    	HostName ige-ssh.u-ga.fr
    	User guilremy
    	IdentityFile ~/.ssh/ige_ssh
    ```
4. Connect to IGE computer and then ADASTRA without needing to login on IGE.

### Important commands

#### Environment variables 
```bash
$HOMEDIR = /lus/home/CT1/c1601279/rguillermin
$WORKDIR = /lus/work/CT1/c1601279/rguillermin
$SCRATCHDIR = /lus/scratch/CT1/c1601279/rguillermin
$STOREDIR = /lus/store/CT1/c1601279/rguillermin
```

#### Current project
To display information about the project (files number, storage, etc):
```bash
myproject -s
```

### Modules
To print activated modules:
```bash
module list
```

### Python environment
To create the Python environment, use the script `env.sh`:
```bash
#!/bin/bash

# Uncomment only if you do NOT source this script.
# set -eu

module purge

module load cpe/24.07
module load cray-python

module list

python3 -m pip install --user --upgrade pip
pip3 install --user --upgrade virtualenv
python3 -m virtualenv ./python_environment
chmod +x ./python_environment/bin/activate
source ./python_environment/bin/activate
python3 -m pip install --upgrade pip
```

Activate the environment:
```bash
source ./python_environment/bin/activate
```

### Install dependencies on ADASTRA

#### Install CARTOPY offline

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

#### Install ESMFpy offline

On ADASTRA, one should manually install `ESMFpy` by downloading and compiling the source code from [earthsystemmodeling.org](https://earthsystemmodeling.org/download/)
```bash
cd esmf-8.8.0/
make info
```

this should output
```
makefile:13: *** Environment variable ESMF_DIR needs to be set to the top ESMF directory. Please see the README file for examples of setting ESMF_DIR.  Stop.
```

One should execute
```bash
echo 'export ESMF_DIR=/lus/home/CT1/c1601279/rguillermin/esmf-8.8.0' >> ~/.bashrc
module purge 
module load cray-libsci_acc/24.07.0
module load intel/2022.1.0
make
```

This will build `ESMF`. Those two modules works together but maybe other combination would work too.

One should check the build with `make all_tests`.

When `EMSF` has been built and checked, one can proceeds to the installation of `esmpy`

```bash
echo 'export ESMFMKFILE=/lus/home/CT1/c1601279/rguillermin/esmf-8.8.0/lib/libO/Linux.gfortran.64.mpiuni.default/esmf.mk' >> ~/.bashrc
cd src/addon/esmpy
python3 -m pip install .
source ~/.bashrc
```

And now one can import `esmpy` with 
```python
import esmpy
```

### Lisa's data
To acces Lisa's data
```bash
cd $STOREDIR
cd ../lweiss/RUN_CROCO/
cd run_swio2_deter_2017_2023_restart/
```

Print the header of a `.nc` file:
```bash
ncdump -h swio_avg.nc 
```

Print its data
```bash
ncview swio_avg.nc
```

### Data transfer
To transfer data from ADASTRA to my computer:
1. compress the files:
    ```bash
    tar -czvf lweiss_scripts.tar.gz ../lweiss/PYTHON/scripts/*/*.py
    ```
2. Disconnect from ADASTRA:
    ```bash
    exit
    ```
3. Copy the data using ige-ssh:
    ```bash
    ssh ige-ssh
    scp rguillermin@adastra.cines.fr:/lus/home/CT1/c1601279/rguillermin/lweiss_scripts.tar.gz .
    exit
    scp -i ~/.ssh/ige_ssh guilremy@ige-ssh.u-ga.fr:lweiss_scripts.tar.gz .
    tar -xzvf lweiss_scripts.tar.gz
    ```

## Satellite data
### Copernicus
We can download datasets from Copernicus by using their python API: [fetch_CMEMS.py](../scripts/copernicus/fetch_CMEMS.py)
