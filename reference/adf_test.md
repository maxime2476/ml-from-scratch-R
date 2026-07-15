# Test de racine unitaire de Dickey-Fuller augmente (ADF)

Regression \\\Delta x_t=\alpha+\beta t+\gamma x\_{t-1}+\sum\_{i=1}^{k}
\delta_i\Delta x\_{t-i}+\varepsilon_t\\ ; la statistique est le t de
\\\gamma\\. \\H_0:\gamma=0\\ (racine unitaire, NON stationnaire).
Valeurs critiques de Dickey-Fuller (non gaussiennes) : rejet si la stat
est TRES negative.

## Usage

``` r
adf_test(x, lags = trunc((length(x) - 1)^(1/3)))
```

## Arguments

- x:

  serie

- lags:

  nombre de retards \\k\\ (defaut : Schwert)

## Value

liste : `statistic`, `lags`.
