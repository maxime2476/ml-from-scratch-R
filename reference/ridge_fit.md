# Régression ridge (forme fermée, éq. 4.2)

Résout \\(X^TX+\lambda I)\beta = X^Ty\\ sur données éventuellement
standardisées, puis re-transforme les coefficients sur l'échelle
d'origine. L'intercept n'est pas pénalisé.

## Usage

``` r
ridge_fit(X, y, lambda, standardize = TRUE, intercept = TRUE)
```

## Arguments

- X:

  matrice de design n x p (sans colonne de constante).

- y:

  réponse.

- lambda:

  pénalité \\\lambda \ge 0\\.

- standardize:

  centrer-réduire les colonnes de X (défaut TRUE).

- intercept:

  inclure un intercept non pénalisé (défaut TRUE).

## Value

liste : `coefficients` (avec intercept si demandé), `beta` (pentes),
`intercept`, `lambda`, `fitted`.
