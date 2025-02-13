# Scripts
Dans ce dossier, je vais regrouper tous les scripts que je vais utiliser à un moment ou à un autre.

## Notes
### [wind_stress.py](wind_stress.py)
- Lignes 47-48: On doit convertir la grille déformée (qui suit les latitudes/longitudes) en une grille géographique pour l'affichage. Pour cela, on va utiliser $u_{geo} = u \cos{\theta} - v \sin{\theta}$ et $v_{geo} = u \sin{\theta} + v \cos{\theta}$ avec $\theta$ l'angle ... . Comme dit dans l'article [Vogt-Vincent et al.](../bibliography/gmd-16-1163-2023.pdf) section 2.2, on a une grille plus large à l'équateur qu'à la limite sud.
- Ligne 56: `projection=ccrs.PlateCarree()` est une projection équivalente à une carte plate où les latitudes et longitudes sont représentées par une grille régulière.
- Lignes 66 et 69: `transform=ccrs.PlateCarree()` permet d'assurer la bonne interprétation des coordonnées.
