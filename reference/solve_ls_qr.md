# Moindres carrés par QR de Householder

Résout \\\min\_\beta \\X\beta - y\\^2\\ via la factorisation QR, en
exploitant la préservation de norme par \\Q\\ : éq. (0.7)-(0.8). La
somme des carrés des résidus est obtenue « gratuitement » comme
\\\\c_2\\^2\\. Voie recommandée quand X est bien conditionnée (préserve
\\\kappa_2(X)\\).

## Usage

``` r
solve_ls_qr(X, y)
```

## Arguments

- X:

  matrice de design n x p, de plein rang colonne.

- y:

  vecteur réponse (longueur n).

## Value

liste : `coefficients`, `fitted`, `residuals`, `rss`, `R` (\\R_1\\),
`Qty` (\\Q^T y\\).
