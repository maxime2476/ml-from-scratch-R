# Quantile conforme (éq. 19.1)

Renvoie le \\\lceil(1-\alpha)(n+1)\rceil\\-ème plus petit score
(\\+\infty\\ si ce rang dépasse n).

## Usage

``` r
conformal_quantile(scores, alpha = 0.1)
```

## Arguments

- scores:

  scores de non-conformité de calibration.

- alpha:

  niveau (couverture visée \\1-\alpha\\).

## Value

le quantile conforme \\\hat q\\.
