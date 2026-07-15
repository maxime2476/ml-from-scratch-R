# Modele de Cox a risques proportionnels (vraisemblance partielle)

Le risque instantane est \\\lambda(t\mid
x)=\lambda_0(t)\exp(x^\top\beta)\\ (proportionnel : le risque de base
\\\lambda_0\\ n'a pas besoin d'etre specifie). On estime \\\beta\\ par
la **vraisemblance partielle** de Cox
\\\prod\_{i:\text{evt}}\exp(x_i^\top\beta)/\sum\_{j\in
R(t_i)}\exp(x_j^\top\beta)\\ (gestion des ex-aequo de Breslow),
maximisee par Newton.

## Usage

``` r
cox_ph(time, event, X)
```

## Arguments

- time, event:

  durees et indicateurs ; @param X matrice de covariables.

## Value

liste : `coefficients`, `se`, `loglik`, `hazard_ratio`.
