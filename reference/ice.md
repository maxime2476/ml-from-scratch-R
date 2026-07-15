# Espérances conditionnelles individuelles (ICE, éq. 15.3)

Le PDP est la moyenne des courbes ICE (colonne `pdp`).

## Usage

``` r
ice(predict_fn, X, feature, grid = NULL, grid_size = 25L)
```

## Arguments

- predict_fn:

  fonction `data.frame -> numeric`.

- X:

  data.frame.

- feature:

  nom de la variable.

- grid:

  grille (sinon régulière).

- grid_size:

  taille par défaut.

## Value

liste : `grid`, `ice` (n x \|grid\|), `pdp` (moyenne des courbes).
