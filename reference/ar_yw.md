# Estimation d'un AR(p) par les equations de Yule-Walker

Resout \\R\phi=r\\ (R matrice de Toeplitz des autocorrelations).
Estimateur de moments, coherent et rapide.

## Usage

``` r
ar_yw(x, order = 1L)
```

## Arguments

- x:

  serie

- order:

  ordre p

## Value

liste : `ar` (coefficients), `var_pred` (variance d'innovation), `mean`.
