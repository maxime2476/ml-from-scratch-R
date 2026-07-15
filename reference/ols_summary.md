# Tableau récapitulatif d'un ajustement MCO

Statistiques t et p-values (éq. 1.8), \\R^2\\ et \\\bar R^2\\ (éq.
1.11-1.12), et test F global (éq. 1.9, cas pente nulle). Reproduit
[`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html).

## Usage

``` r
ols_summary(object)
```

## Arguments

- object:

  objet `ols`.

## Value

liste : `coefficients` (data.frame estimate/se/t/p_value), `r2`,
`adj_r2`, `sigma`, `fstatistic` (value, numdf, dendf, p_value).
