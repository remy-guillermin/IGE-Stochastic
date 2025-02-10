# IGE-Stochastic
Repo pour mon stage à l'IGE sur l'étude des sources d'incertitudes associées à un ensemble de simulation stochastiques.

## Structure du repo
- [bibliography](bibliography/) contient les documents de bilbiographie.
- [book](book/) contient des rapports journaliers et résumés hebdomadaires de mes avancées.

## Connexion à Adastra
Pour me connecter à Adastra je fais 
```bash
ssh adastra
```

Puis je rentre mon mot de passe d'Adastra et je peux activer mon environnement python
```bash
source ./python_environment/bin/activate
```  

## Makefile
- `make` va créer les rendus (slides et rapport).
- `make tree` va créer le fichier `project.tree` représentant l'arborescence du projet.
- `make file.pdf` va compiler le fichier `.tex` associé (report.pdf et slide.pdf disponible).
- `make cleanpdf` va supprimer tous les fichiers pdf compilés.
- `make cleanaux` va supprimer tous les fichiers auxiliaires de compilation.
- `make clean` va executer toutes les commandes `make clean*`.


> [!NOTE]
> À compléter.

## Liens
- [Adastra documentation](https://dci.dci-gitlab.cines.fr/webextranet/)
- [Adastra access documentation](https://dci.dci-gitlab.cines.fr/webextranet/user_support/index.html#adastra-accessing-account-opening)
- [Croco documentation](https://croco-ocean.gitlabpages.inria.fr/croco_doc/)
- 