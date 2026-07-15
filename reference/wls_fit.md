# Moindres carrés pondérés (WLS, éq. 2.5)

Résout \\(X^T W X)^{-1} X^T W y\\ en appliquant l'OLS (QR, Module 0) aux
données transformées \\\sqrt{w_i}\\(x_i, y_i)\\. Reproduit
`lm(weights=)`.

## Usage

``` r
wls_fit(formula, data, weights)
```

## Arguments

- formula:

  formule façon `lm`.

- data:

  data.frame.

- weights:

  vecteur de poids \\w_i = 1/\sigma_i^2\\ (longueur n).

## Value

objet de classe `wls` (mêmes champs qu'`ols`).
