# Pente log-log du RMSE en fonction de n (taux de convergence empirique)

Un estimateur \\\sqrt n\\-consistant a un RMSE \\\propto n^{-1/2}\\,
donc une pente \\\approx -0.5\\ en échelle log-log.

## Usage

``` r
rmse_rate(conv)
```

## Arguments

- conv:

  data.frame issu de `convergence_study`.

## Value

la pente estimée (idéalement ~ -0.5).
