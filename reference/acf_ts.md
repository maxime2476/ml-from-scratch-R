# Fonction d'autocorrelation (ACF)

\\\hat\rho_k = \hat\gamma_k/\hat\gamma_0\\,
\\\hat\gamma_k=\frac1n\sum\_{t} (x_t-\bar x)(x\_{t-k}-\bar x)\\.

## Usage

``` r
acf_ts(x, lag.max = 10L)
```

## Arguments

- x:

  serie

- lag.max:

  retard maximal

## Value

vecteur des autocorrelations (retards 0 a `lag.max`).
