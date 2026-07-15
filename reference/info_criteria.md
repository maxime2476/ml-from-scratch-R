# AIC et BIC d'un modèle ajusté (éq. 6.5, 6.7)

Gère les objets `ols` (Module 1 ; variance comptée comme paramètre,
\\k=p+1\\) et `glm_irls` (Module 3 ; \\k=p\\). Reproduit
[`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
de R.

## Usage

``` r
info_criteria(fit)
```

## Arguments

- fit:

  objet `ols` ou `glm_irls`.

## Value

liste : `aic`, `bic`, `loglik`, `k`, `n`.
