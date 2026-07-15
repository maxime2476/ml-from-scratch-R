# Analyse en composantes principales a noyau (kernel PCA)

PCA dans l'espace de caracteristiques d'un noyau (Module 27) : on centre
la matrice de noyau \\K_c=HKH\\ (\\H=I-\tfrac1n\mathbf 1\mathbf
1^\top\\), on la diagonalise, et l'on projette. Capte des directions
**non lineaires**.

## Usage

``` r
kernel_pca(X, k = 2L, gamma = 1)
```

## Arguments

- X:

  matrice n x p ; @param k nb de composantes ; @param gamma echelle RBF.

## Value

liste : `proj` (n x k), `lambda` (valeurs propres).
