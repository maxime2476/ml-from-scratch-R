# Intervalles de confiance bootstrap

Types : "percentile", "basic" (pivotal), "normal", "bca" (bias-corrected
accelerated, éq. 17.2).

## Usage

``` r
boot_ci(bt, level = 0.95, type = c("percentile", "basic", "normal", "bca"))
```

## Arguments

- bt:

  objet `bootstrap`.

- level:

  niveau de confiance.

- type:

  "percentile", "basic", "normal" ou "bca".

## Value

vecteur (borne inf, borne sup).
