# Regression de Nadaraya-Watson (local constant)

\\\hat m(x_0)=\sum_i K_h(x_0-x_i)\\y_i/\sum_i K_h(x_0-x_i)\\ : moyenne
locale ponderee par le noyau. Biais d'ordre \\h^2\\, mais biais de bord
marque.

## Usage

``` r
nadaraya_watson(x, y, x0, bw)
```

## Arguments

- x, y:

  donnees

- x0:

  points d'evaluation

- bw:

  fenetre

## Value

vecteur des valeurs ajustees en `x0`.
