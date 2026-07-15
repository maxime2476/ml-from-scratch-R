# Trajectoire de la perte d'entraînement/test selon le nombre d'arbres

Trajectoire de la perte d'entraînement/test selon le nombre d'arbres

## Usage

``` r
boost_loss_path(object, data, y)
```

## Arguments

- object:

  objet `boost`.

- data:

  data.frame des prédicteurs.

- y:

  réponse correspondante.

## Value

vecteur de perte (moyenne de \\\ell\\) après m arbres, m = 1..M.
