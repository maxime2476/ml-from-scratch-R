# Test de normalite de Jarque-Bera

\\JB = \tfrac n6\bigl(S^2 + \tfrac14(K-3)^2\bigr)\\, \\S\\ asymetrie,
\\K\\ aplatissement. \\\sim\chi^2_2\\ sous normalite.

## Usage

``` r
jarque_bera(x)
```

## Arguments

- x:

  vecteur (typiquement des residus).

## Value

liste : `statistic`, `skewness`, `kurtosis`, `p_value`.
