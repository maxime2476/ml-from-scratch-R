# Valeurs de Shapley approchées par échantillonnage de permutations

Estimateur de Štrumbelj-Kononenko : pour chaque tirage, une permutation
et une ligne de référence donnent la contribution marginale de chaque
variable.

## Usage

``` r
shapley_permutation(predict_fn, x, X_ref, n_samples = 2000L, seed = NULL)
```

## Arguments

- predict_fn:

  fonction `data.frame -> numeric`.

- x:

  data.frame d'une ligne.

- X_ref:

  data.frame de référence.

- n_samples:

  nombre de tirages.

- seed:

  graine.

## Value

vecteur nommé des valeurs SHAP estimées.
