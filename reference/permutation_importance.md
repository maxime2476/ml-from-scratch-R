# Importance par permutation (éq. 15.6)

Augmentation de la perte quand on permute chaque variable (brisant son
lien avec y). Perte quadratique par défaut.

## Usage

``` r
permutation_importance(
  predict_fn,
  X,
  y,
  loss = function(yh, y) mean((yh - y)^2),
  n_repeat = 10L,
  seed = NULL
)
```

## Arguments

- predict_fn:

  fonction `data.frame -> numeric`.

- X:

  data.frame.

- y:

  réponse.

- loss:

  fonction `(yhat, y) -> perte` (défaut EQM).

- n_repeat:

  répétitions de la permutation (moyennées).

- seed:

  graine.

## Value

vecteur nommé des importances.
