# DBSCAN (clustering base sur la densite)

Un point est **coeur** s'il a au moins `minPts` voisins dans un rayon
`eps`. Les clusters sont des composantes connexes de points-coeur (et
leurs voisins) ; les points isoles sont du **bruit** (label 0). Detecte
des formes arbitraires, sans fixer le nombre de clusters.

## Usage

``` r
dbscan_fit(X, eps, minPts = 5L)
```

## Arguments

- X:

  matrice n x p

- eps:

  rayon

- minPts:

  densite minimale.

## Value

vecteur d'etiquettes (0 = bruit).
