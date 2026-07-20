# Estimation d'un VAR(p) (equation par equation, OLS)

\\Y_t=c+A_1Y\_{t-1}+\dots+A_pY\_{t-p}+\varepsilon_t\\. Chaque equation
est une regression OLS sur la constante et les retards empiles.

## Usage

``` r
var_fit(Y, p = 1L)
```

## Arguments

- Y:

  matrice n x k des series

- p:

  ordre.

## Value

objet `var_scratch` : `B` (coefficients), `Sigma`, `A` (liste des A_l),
etc.
