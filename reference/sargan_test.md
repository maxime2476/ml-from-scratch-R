# Test de suridentification de Sargan

\\n R^2\\ de la regression des residus 2SLS sur TOUS les instruments ;
\\\sim\chi^2\_{q-k}\\ (q instruments, k regresseurs). Rejet =
instruments invalides (correles au terme d'erreur).

## Usage

``` r
sargan_test(y, X, Z)
```

## Arguments

- y:

  reponse

- X:

  design (endogenes inclus)

- Z:

  instruments

## Value

liste : `statistic`, `df`, `p_value`.
