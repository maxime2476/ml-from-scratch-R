# Décomposition QR par réflexions de Householder

Factorise \\X = QR\\ avec \\Q\\ orthogonale (n x n) et \\R\\
triangulaire supérieure (n x p). Implémente les éq. (0.4)-(0.6) de
derivations/00_linalg.qmd. \\Q\\ est accumulée explicitement (\\Q =
H_1\cdots H_p\\) à des fins pédagogiques.

## Usage

``` r
qr_householder(X)
```

## Arguments

- X:

  matrice n x p, \\n \ge p\\.

## Value

liste `Q` (n x n orthogonale), `R` (n x p triangulaire supérieure).
