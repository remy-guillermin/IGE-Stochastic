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
source $WORKDIR/python_environment/bin/activate
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