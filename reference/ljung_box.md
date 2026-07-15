# Test portmanteau de Ljung-Box (autocorrelation residuelle)

\\Q=n(n+2)\sum\_{k=1}^{K}\hat\rho_k^2/(n-k)\sim\chi^2_K\\ sous absence
d'autocorrelation. Utilise sur les residus d'un modele ajuste.

## Usage

``` r
ljung_box(x, lag = 10L, fitdf = 0L)
```

## Arguments

- x:

  serie (ou residus)

- lag:

  nombre de retards K

- fitdf:

  degres de liberte du modele ajuste (a soustraire).

## Value

liste : `statistic`, `df`, `p_value`.
