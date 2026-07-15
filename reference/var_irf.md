# Fonctions de reponse impulsionnelle (IRF) orthogonalisees

Reponse dynamique de chaque variable a un choc unitaire (orthogonalise
par Cholesky de \\\Sigma\\) dans chaque variable, sur `h` periodes.
Coeur de l'analyse macro-econometrique (propagation des chocs).

## Usage

``` r
var_irf(object, h = 10L)
```

## Arguments

- object:

  objet `var_fit` ; @param h horizon.

## Value

tableau (h+1) x k x k : `irf[t, reponse, choc]`.
