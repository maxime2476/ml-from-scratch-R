# Borne de Hoeffding bilatérale (éq. 13.2)

\\\mathbb P(\|\hat R_n - R\|\ge \varepsilon)\le 2e^{-2n\varepsilon^2}\\
pour une perte dans `[0,1]`.

## Usage

``` r
hoeffding_bound(n, eps)
```

## Arguments

- n:

  taille d'échantillon.

- eps:

  écart \\\varepsilon\\.

## Value

la borne supérieure de probabilité.
