# Documentation
## Adastra
### Connexion et configuration
Configuration du fichier `./ssh/config`
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
Après avoir obtenu login sur Adastra
`ssh adastra` sur ma machine, besoin du *password* de l'IGE ainsi que de celui d'Adastra
#### Configuration `ssh` pour accéder aux machines de l'IGE avec une clé
```bash
ssh-keygen -t rsa -b 4096 -C "remy.guillermin1@etu.univ-grenoble-alpes.fr" -f ~/.ssh/ige_ssh
ssh-copy-id -i ~/.ssh/ige_ssh.pub guilremy@ige-ssh.u-ga.fr
```

Après avoir copié la clé publique sur la machine de l'IGE, on peut modifier le fichier `./ssh/config`
```bash
Host ige-ssh
	HostName ige-ssh.u-ga.fr
	User guilremy
	IdentityFile ~/.ssh/ige_ssh
```

Et on peut maintenant se connecter aux machines de l'IGE sans mot de passe et à Adastra avec uniquement le mot de passe d'Adastra

### Commandes importantes
#### Variables d'environnement 
```bash
$HOMEDIR = /lus/home/CT1/c1601279/rguillermin
$WORKDIR = /lus/work/CT1/c1601279/rguillermin
$SCRATCHDIR = /lus/scratch/CT1/c1601279/rguillermin
$STOREDIR = /lus/store/CT1/c1601279/rguillermin
```

#### Projet actuel
Pour afficher les informations du projet (quotas, stockages, etc)
```bash
myproject -s
```

### Modules
Pour afficher les modules actuels
```bash
module list
```

#### Environnement Python
Je peux utiliser la doc d'Adastra et executer le script `env.sh` qui contient
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

Je peux donc activer mon environnement en me connecter et en executant
```bash
source ./python_environment/bin/activate
```


### Données de Lisa
Pour accéder aux données de Lisa, qui sont stockées sur le repertoire de travail, il faut faire
```bash
cd $STOREDIR
cd ../lweiss/RUN_CROCO/
cd run_swio2_deter_2017_2023_restart/
```

Et il est ensuite possible d'afficher le header d'un fichier `.nc` avec
```bash
ncdump -h swio_avg.nc 
```

Et aussi d'afficher les données avec
```bash
ncview swio_avg.nc
```

### Copie de données
Pour transférer les données d'Adastra (ici les scripts de Lisa enregistrés dans [ce dossier](../scripts/lweiss_scripts/)) sur mon ordi je fais dois faire
```bash
(python_environment) [c1601279] rguillermin@login3:~$ tar -czvf lweiss_scripts.tar.gz ../lweiss/PYTHON/scripts/*/*.py
(python_environment) [c1601279] rguillermin@login3:~$ exit
remyguillermin@eduroam-049197 % ssh ige-ssh
guilremy@ige-ssh:~$ scp rguillermin@adastra.cines.fr:/lus/home/CT1/c1601279/rguillermin/lweiss_scripts.tar.gz .
guilremy@ige-ssh:~$ exit
(main) remyguillermin@eduroam-049197 % scp -i ~/.ssh/ige_ssh guilremy@ige-ssh.u-ga.fr:lweiss_scripts.tar.gz .
(main) remyguillermin@eduroam-049197 % tar -xzvf lweiss_scripts.tar.gz
```

## Packages `croco_plot`
Je développe un package python pour me faciliter la tâches afin d'afficher les cartes dont j'aurais besoin.

### Installation
Pour une installation locale, il faut se situer à la racine du projet et exécuter 
```bash
pip install -e modules/
```

### Commande `croco-ipy-load`
Il est possible d'utiliser le script [croco-ipy-load](../scripts/croco-ipy-load.py) afin de directement lancer `iPython` en exécutant la commande `croco-ipy-load`. Pour cela il faut copier le script dans le`/bin/` de l'environnement.
```bash
mkdir -p $VIRTUAL_ENV/bin 
cp scripts/croco-ipy-load.py $VIRTUAL_ENV/bin/croco-ipy-load
chmod +x $VIRTUAL_ENV/bin/croco-ipy-load
source ./python_environment/bin/activate
```