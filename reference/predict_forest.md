# Prédiction d'une forêt / d'un ensemble baggé

Agrège les B arbres : moyenne (régression) ou vote majoritaire
(classification).

## Usage

``` r
predict_forest(object, newdata)
```

## Arguments

- object:

  objet `forest`.

- newdata:

  data.frame des prédicteurs.

## Value

vecteur de prédictions.
