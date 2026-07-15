# Moindres carrés par équations normales et Cholesky

Résout \\X^T X \beta = X^T y\\ en factorisant \\X^T X = L L^T\\ puis par
descente/remontée : éq. (0.10). Rapide mais subit \\\kappa_2(X)^2\\ (cf.
Prop. 0.1) — à réserver aux problèmes bien conditionnés.

## Usage

``` r
solve_ls_chol(X, y)
```

## Arguments

- X:

  matrice de design n x p, de plein rang colonne.

- y:

  vecteur réponse (longueur n).

## Value

liste : `coefficients`, `fitted`, `residuals`, `rss`, `L`.
