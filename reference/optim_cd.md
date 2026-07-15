# Coordinate descent générique

Minimise cycliquement selon chaque coordonnée (éq. 0.19), en supposant
fourni le minimiseur 1-D exact `argmin_coord(x, j)` (structure «
sous-problème 1-D fermé », qui deviendra le soft-thresholding au Module
4).

## Usage

``` r
optim_cd(argmin_coord, x0, max_sweep = 1000L, tol = 1e-09, f = NULL)
```

## Arguments

- argmin_coord:

  fonction : `argmin_coord(x, j)` renvoie la valeur optimale de la
  coordonnée j, les autres coordonnées de x étant fixées.

- x0:

  point initial.

- max_sweep:

  nombre maximal de balayages complets.

- tol:

  seuil d'arrêt sur la variation d'un balayage.

- f:

  (optionnel) fonction objectif.

## Value

liste : `par`, `sweeps`, `value`.
