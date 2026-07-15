# Prédiction d'un MLP

Prédiction d'un MLP

## Usage

``` r
predict_mlp(model, newdata, type = c("response", "class"))
```

## Arguments

- model:

  objet `mlp`.

- newdata:

  matrice des prédicteurs.

- type:

  "response" (valeur/probabilité) ou "class" (0/1, log-loss).

## Value

vecteur/matrice de prédictions.
