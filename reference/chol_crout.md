# Factorisation de Cholesky (algorithme de Crout)

Factorise une matrice SPD \\A = L L^T\\, \\L\\ triangulaire inférieure à
diagonale positive. Implémente les formules (0.9) de
derivations/00_linalg.qmd.

## Usage

``` r
chol_crout(A)
```

## Arguments

- A:

  matrice carrée symétrique définie positive.

## Value

le facteur `L` (triangulaire inférieur).
