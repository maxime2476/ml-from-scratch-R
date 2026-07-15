# Dépendance partielle (PDP, éq. 15.2)

Dépendance partielle (PDP, éq. 15.2)

## Usage

``` r
pdp(predict_fn, X, feature, grid = NULL, grid_size = 25L)
```

## Arguments

- predict_fn:

  fonction `data.frame -> numeric`.

- X:

  data.frame des données.

- feature:

  nom de la variable d'intérêt.

- grid:

  grille de valeurs (sinon régulière sur l'étendue observée).

- grid_size:

  taille de grille par défaut.

## Value

data.frame `grid`, `pdp`.
