# Convolution 2D (correlation croisee, mode "valid")

Pour chaque filtre \\f\\,
\\Y\_{ij}^f=b_f+\sum\_{a,b}X\_{i+a-1,j+b-1}K\_{ab}^f\\. Partage de poids
: le meme noyau \\K^f\\ est applique en toute position.

## Usage

``` r
conv2d(X, K, b = rep(0, dim(K)[3]))
```

## Arguments

- X:

  image d'entree (matrice H x W).

- K:

  noyaux, tableau (kh x kw x F) ; @param b biais (longueur F).

## Value

liste : `out` (Hout x Wout x F), `cache`.
