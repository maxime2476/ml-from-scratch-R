# Moindres carrés de norme minimale par SVD

Renvoie \\\hat\beta\_{\min} = X^+ y\\ (éq. 0.13), qui minimise
\\\\X\beta - y\\\\ et, en cas de rang déficient, a la plus petite norme
parmi tous les minimiseurs (Prop. 0.5). Solution propre à la colinéarité
parfaite et au cas \\p \> n\\.

## Usage

``` r
solve_ls_svd(X, y)
```

## Arguments

- X:

  matrice de design n x p (rang quelconque).

- y:

  vecteur réponse (longueur n).

## Value

liste : `coefficients`, `fitted`, `residuals`, `rss`, `rank`, `kappa`.
