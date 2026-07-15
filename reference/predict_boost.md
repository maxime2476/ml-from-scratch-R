# Prédiction d'un modèle de boosting

Prédiction d'un modèle de boosting

## Usage

``` r
predict_boost(
  object,
  newdata,
  type = c("response", "link", "class"),
  n_trees = object$M
)
```

## Arguments

- object:

  objet `boost`.

- newdata:

  data.frame des prédicteurs.

- type:

  "response" (valeur / probabilité), "link" (score F) ou "class" (0/1,
  log-loss uniquement).

- n_trees:

  nombre d'arbres à utiliser (défaut : tous).

## Value

vecteur de prédictions.
