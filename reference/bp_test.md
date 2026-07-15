# Test de Breusch-Pagan (heteroscedasticite)

Version studentisee de Koenker : \\n R^2\\ de la regression des residus
au carre sur les regresseurs. \\H_0\\ : homoscedasticite. Statistique
\\\sim \chi^2\_{p-1}\\.

## Usage

``` r
bp_test(formula, data)
```

## Arguments

- formula:

  formule du modele

- data:

  data.frame

## Value

liste : `statistic`, `df`, `p_value`.
