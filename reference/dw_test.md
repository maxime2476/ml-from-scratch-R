# Statistique de Durbin-Watson (autocorrelation d'ordre 1)

\\DW = \sum\_{t\ge2}(e_t-e\_{t-1})^2/\sum_t e_t^2 \approx
2(1-\hat\rho)\\. Proche de 2 : pas d'autocorrelation ; proche de 0 :
positive.

## Usage

``` r
dw_test(formula, data)
```

## Arguments

- formula, data:

  cf. `bp_test`.

## Value

liste : `statistic`, `rho` (autocorrelation implicite).
