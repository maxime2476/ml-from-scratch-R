# Résolution d'un système triangulaire inférieur par descente

Résout \\L x = b\\ avec \\L\\ triangulaire inférieure. Brique de la
résolution des équations normales par Cholesky, éq. (0.10) de
derivations/00_linalg.qmd.

## Usage

``` r
forward_substitution(L, b)
```

## Arguments

- L:

  matrice triangulaire inférieure (n x n), diagonale non nulle.

- b:

  vecteur second membre (longueur n).

## Value

le vecteur solution x.
