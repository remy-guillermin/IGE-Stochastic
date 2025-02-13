# Documentation

## Adastra

### Connexion et configuration

#### Configuration du fichier `./ssh/config`
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

Après avoir obtenu le login sur Adastra, utilisez la commande suivante pour vous connecter :
```bash
ssh adastra
```
Vous aurez besoin du mot de passe de l'IGE ainsi que de celui d'Adastra.

#### Configuration `ssh` pour accéder aux machines de l'IGE avec une clé
1. Générer une clé SSH :
    ```bash
    ssh-keygen -t rsa -b 4096 -C "remy.guillermin1@etu.univ-grenoble-alpes.fr" -f ~/.ssh/ige_ssh
    ```
2. Copier la clé publique sur la machine de l'IGE :
    ```bash
    ssh-copy-id -i ~/.ssh/ige_ssh.pub guilremy@ige-ssh.u-ga.fr
    ```
3. Modifier le fichier `./ssh/config` :
    ```bash
    Host ige-ssh
    	HostName ige-ssh.u-ga.fr
    	User guilremy
    	IdentityFile ~/.ssh/ige_ssh
    ```
4. Se connecter aux machines de l'IGE sans mot de passe et à Adastra avec uniquement le mot de passe d'Adastra.

### Commandes importantes

#### Variables d'environnement 
```bash
$HOMEDIR = /lus/home/CT1/c1601279/rguillermin
$WORKDIR = /lus/work/CT1/c1601279/rguillermin
$SCRATCHDIR = /lus/scratch/CT1/c1601279/rguillermin
$STOREDIR = /lus/store/CT1/c1601279/rguillermin
```

#### Projet actuel
Pour afficher les informations du projet (quotas, stockages, etc) :
```bash
myproject -s
```

### Modules
Pour afficher les modules actuels :
```bash
module list
```

### Environnement Python
Pour configurer l'environnement Python, utilisez le script `env.sh` :
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

Pour activer l'environnement :
```bash
source ./python_environment/bin/activate
```

### Données de Lisa
Pour accéder aux données de Lisa :
```bash
cd $STOREDIR
cd ../lweiss/RUN_CROCO/
cd run_swio2_deter_2017_2023_restart/
```

Pour afficher le header d'un fichier `.nc` :
```bash
ncdump -h swio_avg.nc 
```

Pour afficher les données :
```bash
ncview swio_avg.nc
```

### Copie de données
Pour transférer les données d'Adastra sur votre ordinateur :
1. Compresser les scripts de Lisa :
    ```bash
    tar -czvf lweiss_scripts.tar.gz ../lweiss/PYTHON/scripts/*/*.py
    ```
2. Se déconnecter d'Adastra :
    ```bash
    exit
    ```
3. Copier les données sur votre machine via ige-ssh :
    ```bash
    ssh ige-ssh
    scp rguillermin@adastra.cines.fr:/lus/home/CT1/c1601279/rguillermin/lweiss_scripts.tar.gz .
    exit
    scp -i ~/.ssh/ige_ssh guilremy@ige-ssh.u-ga.fr:lweiss_scripts.tar.gz .
    tar -xzvf lweiss_scripts.tar.gz
    ```

## Packages `croco_plot`
Je développe un package python pour faciliter l'affichage des cartes.

### Installation
Pour une installation locale, exécutez :
```bash
pip install -e modules/
```

### Commande `croco-ipy-load`
Pour utiliser le script `croco-ipy-load` :
1. Copier le script dans le `/bin/` de l'environnement :
    ```bash
    mkdir -p $VIRTUAL_ENV/bin 
    cp scripts/croco-ipy-load.py $VIRTUAL_ENV/bin/croco-ipy-load
    chmod +x $VIRTUAL_ENV/bin/croco-ipy-load
    source ./python_environment/bin/activate
    ```

#### Arguments
- `--clear` : exécute `clear` avant de lancer l'instance `iPython`