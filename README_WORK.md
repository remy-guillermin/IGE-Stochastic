# $WORKDIR

Ici on a l'ensemble de mes données après post traitement (retrait des NaN, calcul de la MLD/tension du vent).

## CROCO

Ce dossier contient le code source `croco` qui a été utilisé juste pour me familiariser avec le modèle.

## CROCO_dev_sto_CD

Ce dossier contient une copie du code `croco` qui a été modifié par Lisa afin d'implémenter une perturbation sur la tension du vent. Je m'en suis servi afin de développer la perturbation du mélange vertical.

## CROCO_dev_sto_GLS

Ce dossier contient le code `croco` avec l'implémentation que j'ai fait pour la perturbation du mélange vertical.

## DATASETS

Ce dossier contient les données des simulations après post traitement et contient uniquement des séries temporelles des valeurs moyennées dans chaque zone pour les traceurs (Salinité, Température et Surface Libre)

Il y a trois dossiers contenant les valeurs de la simulation, de la réanalyse GLORYS et de différentes sourcesz d'observation.

## grid 

Ce dossier contient simplement les différents fichiers de grille. 

## MLD

Ce dossier contient les fichiers `.nc` de données post traitées pour le calcul de la profondeur de la couche de mélange. Il y a aussi bien la valeur interpolée linéairement que la valeur un niveau au dessus.

## MLD_old

Ce dossier contient la première itération de fichiers `.nc` de données post traitées pour le calcul de la profondeur de la couche de mélange. Les valeurs sont faussées à cause de la méthode de calcul qui ne prenaient pas en compte des variations locales de densité à l'origine.

## NaN_CORRECTED


## NaN_MERGED

Même principe que précédemment sauf qu'ici il y a un fichier par membre qui contient les trois ans. On y retrouve aussi les fichiers de moyenne et d'écart-type d'ensemble.

## OBS

Ce dossier contient les observations venant des différents produits du CMEMS traitées pour mes besoins (retrait des variables inutiles ainsi que des zones non utilisées)

## RAW

Ce dossier contient les données d'observations brutes.

## REGRIDDED

Ce dossier contient les observations interpolées sur la même grille que celle de nos simulations.

# WINDSTR

Ce dossier contient les fichiers `.nc` de données post traitées pour le calcul de la tension du vent. Il y a aussi bien la valeur selon la direction u, que v ainsi que la "norme".