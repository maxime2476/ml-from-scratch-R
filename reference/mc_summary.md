# Résumé d'une étude de simulation : biais, RMSE, variance, avec erreurs MC

Résumé d'une étude de simulation : biais, RMSE, variance, avec erreurs
MC

## Usage

``` r
mc_summary(estimates, truth)
```

## Arguments

- estimates:

  vecteur des R estimations \\\hat\theta\\.

- truth:

  valeur vraie \\\theta\\.

## Value

liste : `R`, `mean`, `bias` (+ `bias_se`), `rmse` (+ `rmse_se`),
`variance`, `empirical_se`.
