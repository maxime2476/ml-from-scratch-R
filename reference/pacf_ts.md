# Fonction d'autocorrelation partielle (PACF, recursion de Durbin-Levinson)

Fonction d'autocorrelation partielle (PACF, recursion de
Durbin-Levinson)

## Usage

``` r
pacf_ts(x, lag.max = 10L)
```

## Arguments

- x:

  serie

- lag.max:

  retard maximal

## Value

vecteur des autocorrelations partielles (retards 1 a `lag.max`).
