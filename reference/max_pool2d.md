# Max-pooling 2D (sous-echantillonnage par le maximum)

Reduit chaque bloc `pool` x `pool` a son maximum (invariance locale a la
translation, reduction de dimension).

## Usage

``` r
max_pool2d(X, pool = 2L)
```

## Arguments

- X:

  carte de caracteristiques (H x W) ; @param pool taille du bloc.

## Value

liste : `out`, `cache` (positions des maxima).
