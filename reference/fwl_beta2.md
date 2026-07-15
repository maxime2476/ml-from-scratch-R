# Coefficient FWL : régression de M1 y sur M1 X2 (éq. 1.10)

Renvoie \\\hat\beta_2\\ obtenu en résidualisant y et X2 par rapport à X1
(théorème de Frisch-Waugh-Lovell). Sert à la vérification numérique : ce
vecteur doit égaler le sous-vecteur \\\hat\beta_2\\ de la régression
complète.

## Usage

``` r
fwl_beta2(y, X1, X2)
```

## Arguments

- y:

  réponse (vecteur).

- X1:

  bloc à partialiser (matrice n x p1).

- X2:

  bloc d'intérêt (matrice n x p2).

## Value

le vecteur \\\hat\beta_2\\.
