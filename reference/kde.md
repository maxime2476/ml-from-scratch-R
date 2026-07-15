# Estimation de densite par noyau (KDE)

\\\hat f(x)=\frac1{nh}\sum_i K\\\bigl(\frac{x-x_i}{h}\bigr)\\, noyau
gaussien.

## Usage

``` r
kde(x, grid, bw)
```

## Arguments

- x:

  donnees

- grid:

  points d'evaluation

- bw:

  fenetre h

## Value

vecteur des densites estimees sur `grid`.
