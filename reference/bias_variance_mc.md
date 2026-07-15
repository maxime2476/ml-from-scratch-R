# Estimation Monte Carlo de la décomposition biais-variance (éq. 6.1)

Sur un DGP connu, estime irréductible / biais² / variance de la
prédiction en un point de test x0, pour un ajusteur donné.

## Usage

``` r
bias_variance_mc(gen_y, fit_pred, f0, sigma2, R = 5000L)
```

## Arguments

- gen_y:

  fonction `() -> y` générant une réponse (design X0 fixe, bruit
  aléatoire).

- fit_pred:

  fonction `(y) -> prédiction en x0` (scalaire).

- f0:

  vraie valeur \\f(x_0)\\.

- sigma2:

  variance du bruit (partie irréductible).

- R:

  nombre de réplications.

## Value

liste : `irreducible`, `bias2`, `variance`, `mse_pred`, `total`.
