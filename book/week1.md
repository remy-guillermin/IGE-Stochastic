# Semaine du Lundi 3 Fevrier

Première semaine de stage ! 

## Lundi
Aujourd'hui c'était journée administrative: Le matin réunion nouveaux arrivants et l'après midi RDV avec l'administration pour récupérer les clés. 

Lisa et Jean-Michel m'ont donné 3 articles comme base de travail: 
- [Article sur les simulations stochastiques](../bibliography/Brankart%20et%20al.%20-%202015%20-%20A%20generic%20approach%20to%20explicit%20simulation%20of%20uncer.pdf): *Une description de la méthode probabiliste et ensembliste sur laquelle se base la méthodologie de ton stage. On souhaite appliquer cette méthode à notre modèle régional afin de quantifier les incertitudes associées à la circulation de surface du sud-ouest de l'Océan Indien.*
- [Article sur les processus océaniques Indien](../bibliography/gmd-16-1163-2023.pdf) *Une review décrivant le contexte physique des processus océaniques de l'Océan Indien. A partir de cette review, il serait intéressant de dégager les processus affectant la circulation dans la partie Sud-Ouest de l'océan Indien et les questions scientifiques qui en découlent (étendue de notre zone d'étude en pj).*
- [Article sur une configuration similaire de Croco](../bibliography/os-17-1677-2021.pdf) *Une description d'une configuration régionale du modèle d'océan CROCO sur le Sud-Ouest de l'océan Indien, qui est construite d'une façon analogue à celle que tu utiliseras dans ton stage, associée à une analyse climatologique des sorties de modèle.*

On a aussi fait la demande d'accès à [Adastra](https://dci.dci-gitlab.cines.fr/webextranet/)

## Mardi
Aujourd'hui pas grand chose, j'ai commencé à travailler en profondeur l'[article sur les simulations stochastiques](../bibliography/Brankart%20et%20al.%20-%202015%20-%20A%20generic%20approach%20to%20explicit%20simulation%20of%20uncer.pdf).

> [!IMPORTANT]
> NEMO a une résolution verticale fixe et est donc moins précis sur les zones a faible bathymétrie (relief marin).
> 
> CROCO lui a une résolution verticale qui suit la bathymétrie ce qui permet d'être mieux résolu sur les côtes notamment.

## Mercredi
### Questions sur l'article de Jean-Michel
> **À quoi sert la *correlation timescale* ?**
>
> Il s'agit de l'échelle de temps de corrélation, par exemple une Gaussienne plutôt plate qui décrit combien de temps une mesure en un point va être correlée à un autre mesure en ce même point à un temps différents.


> **Pourquoi résoudre une équation elliptique pour faire un filtre spatial ?**
> 
> Un équation elliptique est une PDE de la forme $L^2 \Delta \psi = w$ avec $w$ un bruit blanc, donc non-corrélé spatiallement. Quand on résoud cette équation on va obtenir une fonction $\psi$ qui sera corrélée dans l'espace.


> **Pourquoi est-on intéressé par l'évolution de la biodiversité marine dans les modèles numériques physiques ?**
>
> L'importance n'est pas les incertitudes causées par la biodiversité marine mais les incertitudes **sur** la biodiversité marine. Par exemple dans l'article de Jean-Michel uniquement les deux plus petites échelles de biodiversité sont résolues.

## Jeudi
Aujourd'hui on a finalisé la configuration de l'accès à Adastra, je peux donc me connecter à Adastra et utiliser Python.