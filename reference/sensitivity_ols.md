# Analyse de sensibilité complète d'un effet OLS

Analyse de sensibilité complète d'un effet OLS

## Usage

``` r
sensitivity_ols(fit, treatment, q = 1, alpha = 0.05)
```

## Arguments

- fit:

  objet `ols` (Module 1).

- treatment:

  nom du régresseur d'intérêt.

- q:

  fraction de réduction pour la robustness value.

- alpha:

  seuil pour la robustness value significative.

## Value

liste : `estimate`, `se`, `t`, `df`, `r2yd`, `rv_q`, `rv_qa`.
