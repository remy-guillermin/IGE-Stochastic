# Semaine du Lundi 10 Février

## Lundi
Adastra est en panne aujourd'hui, donc je fais de la lecture. J'ai lu et annoté [l'article sur l'analyse climatologique](../bibliography/gmd-16-1163-2023.pdf) et je vais commencer [celui sur le contexte physique](../bibliography/os-17-1677-2021.pdf).

## Mardi
Adastra est toujours en panne ce matin, donc je continue de travailler sur le contexte dynamique de la région.

- **Connexion entre l'océan Pacifique et l'océan Indien:** Le flot indonésien traversant (*Indonesian ThroughFlow* ITF) injecte de l'eau douce et chaude dans l’océan Indien tropical du sud. Ces eaux, plus légères que celles plus au sud, créent un gradient nord-sud de densité et de pression, alimentant différents courants.
  - *Le courant géostrophique de surface*, large et dirigé vers l’est entre 16 et 32° S, entre Madagascar et l'Australie.
  - *Le courant de Leeuwin*, dirigé vers l'Antarctique et longeant l'Australie.
- **Courant d'Agulhas:** Ce courant de surface et d'eau moyennement profonde longe la pointe de l'Afrique, déplaçant rapidement de la chaleur vers le sud. Il contribue à 30% de l'exportation de chaleur de l'océan Indien.
- **Le courant de Leeuwin:** Moins contributif à l'export de chaleur que le courant d'Agulhas, il injecte une grande variété de tourbillons méso-échelle transportant chaleur et quantité de mouvement.
- **Les cellules:** Deux cellules connectent différentes zones de remontées d'eau dans le sud et le nord de l'océan Indien, jouant un rôle majeur dans la régulation de l'équilibre moyen, tropical et interannuel de chaleur.
  - *La cellule sub-tropicale* emmène de l'eau froide jusqu'au nord de l'océan Indien.
  - *La cellule cross-équatoriale* ramène l'eau froide du pôle jusqu'au milieu de l'océan Indien.
- **Les eaux modales[^1]:** 
  - *À des profondeurs intermédiaires (500-2000 m)*, les eaux modales de l'océan Austral entrent dans l'océan Indien. En remontant, elles se mélangent avec les eaux de surface plus légères et entament un mouvement de remontée vers la surface dans différentes régions au nord de 10° S, puis sont ramenées vers le sud par le transport d'Ekman[^2].
  - *En profondeur (2000-4000 m)*, les eaux modales se mélangent avec les eaux plus denses et rejoignent le courant d'eau profonde allant vers le sud, contribuant aux eaux profondes d'Antarctique se déplaçant vers le nord.
- **L'écoulement cross-équatorial:** Accompli aux profondeurs abyssales, ainsi qu'avec le Courant de Côte d'Afrique de l'Est (EACC), le Courant Somalien inversant saisonnièrement et le transport d'Ekman vers le sud.
- **El Niño:** Anomalie positive et pseudo-périodique (2 à 7 ans) des températures de l'océan Pacifique, entraînant un vent anticyclonique dans l'océan Indien du Sud-Est.
  - *El Niño Southern Oscillation (ENSO)* relie El Niño et l'oscillation australe de la pression atmosphérique, alimentant les variations de SSS (*Sea Surface Salinity*) et maintenant le réchauffement de l'océan Indien tropical par la propagation vers l'est d'ondes de Rossby.
- **La Niña:** Anomalie négative pseudo-périodique, non corrélée à El Niño, des températures de l'océan Pacifique, entraînant un vent cyclonique.
  - ENSO prolonge les anomalies de température dans l'océan Indien, de manière similaire à El Niño.

[^1]: Masses d’eau homogènes formées par le mélange et la convection dans certaines régions des océans, souvent associées à des gyres subtropicaux.
[^2]: Mouvement des eaux de surface dû à l'action des contraintes des vents et de la force de Coriolis.

![](indian_ocean_circulation.png)
*Figure 1: Vue schématique des phénomènes clés dans l'Océan Indien. https://doi.org/10.5194/os-17-1677-2021*

J'ai commencé à tracer différentes cartes, pour le moment des cartes de vent et de bathymétrie. Demain, je vais essayer de faire une carte d'EKE et travailler sur les zones de Lisa qui sont dans [ce script](../scripts/lweiss_mod/KE/eke_avg_croco.py).

## Mercredi
J'ai bien avancé sur les figures aujourd'hui. J'ai des scripts pour afficher les vitesses, vorticité, hélicité ($u \times \omega$) ainsi que les SST, SSH et SSS. Mon package est à jour par rapport à mes figures.

> [!IMPORTANT]
> Demain, il faut modifier la `cmap` de l'EKE car elle est en log mais pas les ticks.

## Jeudi
Aujourd'hui, je commence l'analyse des zones, avec en premier l'EKE de surface et la SLA. J'ai fait des moyennes glissantes pour la SLA en plus.

### Résumé des commits d'aujourd'hui
- Ajouté des clarifications et amélioré la structure des compte rendus des semaines précédentes.
- Mis à jour les figures et les scripts pour afficher les vitesses, vorticité, hélicité, SST, SSH et SSS.
- Travaillé sur les cartes de vent et de bathymétrie.
- Préparé une carte d'EKE et commencé à analyser les zones avec l'EKE de surface et la SLA.

## Vendredi 
Aujourd'hui, j'ai continué à travailler sur mes scripts pour les séries temporelles, l'énergie et les anomalies de la mer.

- **Séries temporelles:** J'ai amélioré les scripts pour générer des séries temporelles plus précises des différentes variables océanographiques.
- **Énergie:** J'ai analysé les flux d'énergie dans l'océan Indien, en me concentrant sur les échanges entre les différentes couches de l'océan.
- **Anomalies de la mer:** J'ai travaillé sur la détection et l'analyse des anomalies de la mer, en particulier les anomalies de température et de salinité de surface.
- **Cartes de vent et de bathymétrie:** J'ai finalisé les cartes de vent et de bathymétrie et commencé à travailler sur les cartes d'EKE.
- **Analyse des zones:** J'ai approfondi l'analyse des zones avec l'EKE de surface et la SLA, en utilisant des moyennes glissantes pour la SLA.

### Résumé des commits d'aujourd'hui
- Amélioration des scripts pour les séries temporelles.
- Analyse des flux d'énergie dans l'océan Indien.
- Détection et analyse des anomalies de la mer.
- Finalisation des cartes de vent et de bathymétrie.
- Travail sur les cartes d'EKE et analyse des zones avec l'EKE de surface et la SLA.

## TODO pour Lundi
- Réfléchir à un plan pour le premier chapitre de mon rapport.