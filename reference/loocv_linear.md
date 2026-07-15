# LOOCV fermé pour la régression linéaire (éq. 6.3)

Calcule \\\mathrm{CV}\_n = \frac1n\sum_i
(\hat\varepsilon_i/(1-h\_{ii}))^2\\ en UN seul ajustement, via les
leviers \\h\_{ii}\\ (Théorème 6.2).

## Usage

``` r
loocv_linear(X, y)
```

## Arguments

- X:

  matrice de design (constante incluse).

- y:

  réponse.

## Value

liste : `cv` (erreur LOOCV moyenne), `loo_resid`, `h` (leviers).
