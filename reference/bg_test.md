# Test de Breusch-Godfrey (autocorrelation d'ordre p, LM)

Regression des residus sur les regresseurs ET \\p\\ retards des residus
; \\n R^2 \sim \chi^2_p\\. Plus general que Durbin-Watson (ordre eleve,
regresseurs retardes admis).

## Usage

``` r
bg_test(formula, data, order = 1)
```

## Arguments

- formula, data:

  cf. `bp_test`

- order:

  ordre du retard

## Value

liste : `statistic`, `df`, `p_value`.
