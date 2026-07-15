# Distance de l'origine à l'enveloppe convexe de colonnes (QP exact)

Résout \\\min\_{\lambda\ge 0,\\ \mathbf 1^\top\lambda=1}\\V\lambda\\\\
par programmation quadratique. Sert au test de séparabilité (théorème de
Gordan : un étiquetage est linéairement séparable ssi 0 n'est PAS dans
l'enveloppe convexe des \\y_i\\\tilde x_i\\).

## Usage

``` r
.min_dist_hull(V)
```

## Arguments

- V:

  matrice (colonnes = points).

## Value

la distance minimale de l'origine à l'enveloppe convexe.
