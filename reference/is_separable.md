# Séparabilité linéaire d'un étiquetage (théorème de Gordan)

Séparabilité linéaire d'un étiquetage (théorème de Gordan)

## Usage

``` r
is_separable(X, y, tol = 1e-06)
```

## Arguments

- X:

  matrice n x d des points.

- y:

  étiquettes dans {-1, +1} (longueur n).

- tol:

  seuil de distance pour déclarer 0 hors de l'enveloppe.

## Value

TRUE si {(x_i, y_i)} est linéairement séparable (avec biais).
