# Jackknife (delete-one) — le jackknife infinitésimal EST la fonction d'influence

Pour un estimateur lisse, \\\hat\theta\_{(-i)}-\hat\theta \approx
-\mathrm{IC}\_i/n\\, d'où \\\widehat{\operatorname{Var}}\_{\text{jack}}
= \frac{n-1}{n}\sum_i (\hat\theta\_{(-i)}-\bar\theta\_{(\cdot)})^2
\approx \frac1{n^2}\sum \mathrm{IC}\_i^2\\.

## Usage

``` r
jackknife(data, stat_fn)
```

## Arguments

- data:

  data.frame (une ligne = une observation).

- stat_fn:

  fonction `data -> theta` (scalaire).

## Value

liste : `estimate`, `bias`, `var`, `values`.
