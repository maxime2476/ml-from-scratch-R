# Moindres carrés généralisés (GLS, éq. 2.4) avec Omega connue

Résout \\(X^T \Omega^{-1} X)^{-1} X^T \Omega^{-1} y\\ en appliquant
l'OLS aux données transformées \\P^{-1}(X, y)\\ où \\\Omega = P P^T\\
(Cholesky, Module 0). Théorème d'Aitken (Th. 2.3).

## Usage

``` r
gls_fit(formula, data, Omega)
```

## Arguments

- formula:

  formule façon `lm`.

- data:

  data.frame.

- Omega:

  matrice n x n SPD (structure de covariance des erreurs, à un facteur
  d'échelle sigma^2 près).

## Value

objet de classe `gls_fit` (mêmes champs qu'`ols`).
