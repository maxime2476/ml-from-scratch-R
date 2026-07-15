# Résolution d'un système triangulaire supérieur par remontée

Résout \\U x = b\\ avec \\U\\ triangulaire supérieure. Brique de la
résolution MCO par QR (remontée sur \\R_1\\, éq. 0.8) et par Cholesky
(éq. 0.10).

## Usage

``` r
back_substitution(U, b)
```

## Arguments

- U:

  matrice triangulaire supérieure (n x n), diagonale non nulle.

- b:

  vecteur second membre (longueur n).

## Value

le vecteur solution x.
