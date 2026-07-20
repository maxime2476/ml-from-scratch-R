# Analyse en composantes independantes (FastICA)

Separe un melange \\X=SA^\top\\ en **sources independantes** \\S\\, en
maximisant la **non-gaussianite** (negentropie, contraste \\g=\tanh\\)
par iteration de point fixe, apres **blanchiment**. Va au-dela de la PCA
(qui ne fait que decorreler) : resout le "cocktail party".

## Usage

``` r
ica_fastica(X, n_comp = ncol(X), iter = 300L, tol = 1e-09)
```

## Arguments

- X:

  melange (n x d)

- n_comp:

  nb de composantes

- iter, tol:

  arret.

## Value

liste : `S` (sources estimees), `W` (matrice de separation).
