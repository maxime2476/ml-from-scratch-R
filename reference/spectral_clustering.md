# Clustering spectral

Construit un graphe de similarite (noyau RBF), calcule le **Laplacien
normalise** \\L=I-D^{-1/2}WD^{-1/2}\\, plonge les points dans les `k`
premiers vecteurs propres, puis applique k-means. Separe des clusters
**non convexes** (spirales, cercles) que le k-means echoue a distinguer.

## Usage

``` r
spectral_clustering(X, k, gamma = 1)
```

## Arguments

- X:

  matrice n x p

- k:

  nombre de clusters

- gamma:

  echelle RBF.

## Value

vecteur d'assignation.
